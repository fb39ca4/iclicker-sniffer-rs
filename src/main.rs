#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(unused_mut)]
#![allow(dead_code)]

extern crate core;
extern crate num;
extern crate rustfft;
extern crate itertools;

use std::error::Error;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::process::{Command, Stdio};
use std::option;
use std::vec::Vec;
use std::boxed::Box;
use std::fmt::Display;

use itertools::Itertools;

use num::complex::Complex32;
use num::complex::Complex;
use num::traits::Float;
use num::traits::Zero;
use core::ops::Neg;

use std::path::Path;
use std::fs::File;

struct SquelchedSamples<I> {
    count: u8,
    iterator: I,//Iterator<Item=std::io::Result<u8>>,
}

impl<I> SquelchedSamples<I> where I: Iterator<Item=Complex32> {
    fn treshold(&self, x: Complex32) -> bool {
        x.norm() > (0.1 as f32)
    }
}

impl<I> Iterator for SquelchedSamples<I> where I: Iterator<Item=Complex32> {
    type Item = Complex32;
    fn next(&mut self) -> Option<Complex32> {
        loop {
            //return self.iterator.next();
            let opt = self.iterator.next();
            match opt {
                Some(value) => (),
                None => {
                    return None;
                }
            }
            let x = opt.unwrap();
            if self.treshold(x) {
                self.count = 32;
            }
            else {
                if self.count > 0 {
                    self.count -= 1;
                }
            }
            if self.count > 0 {
                return Some(x);
            }
        }
    }
}

struct QuadratureDemod<T, I> where T: Neg<Output=T> + Copy + Float + Default, I: Iterator<Item=Complex<T>> {
    iterator: I,
    previous: Complex<T>,
    file: Option<Box<Write>>,
}

impl<T, I> QuadratureDemod<T, I> where T: Neg<Output=T> + Copy + Float + Default, I: Iterator<Item=Complex<T>> {
    fn new(samples: I, file: Option<Box<Write>>) -> QuadratureDemod<T, I> {
        QuadratureDemod {
            iterator: samples,
            previous: Complex::default(),
            file: file,
        }
    }
}

impl<T, I> Iterator for QuadratureDemod<T, I> where T: Neg<Output=T> + Copy + Float + Default + Display, I: Iterator<Item=Complex<T>> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        let new_sample = match self.iterator.next() {
            Some(value) => value,
            None => {
                return None;
            }
        };
        let conj_prev = self.previous.conj();
        self.previous = new_sample;
        let freq = (new_sample * conj_prev).to_polar().1;
        match self.file {
            Some(ref mut file) => {
                write!(file, "{}\n", freq).unwrap();
            },
            None => (),
        }
        Some(freq)
    }
}

struct FFTDemod<I> {
    window_size: usize,
    iterator: I,
    window: std::collections::VecDeque<Complex32>,
    fft_input: Vec<Complex32>,
    fft_output: Vec<Complex32>,
    fft: rustfft::FFT<f32>,
    file: Option<Box<Write>>,
}

impl<I> FFTDemod<I> where I: Iterator<Item=Complex32> {
    fn new(samples: I, window_size: usize, file: Option<Box<Write>>) -> FFTDemod<I> {
        let mut demod_struct = FFTDemod {
            iterator: samples,
            window_size: window_size,
            window: std::collections::VecDeque::with_capacity(window_size),
            fft_input: vec![Complex32::default(); window_size],
            fft_output: vec![Complex32::default(); window_size],
            fft: rustfft::FFT::new(window_size, false),
            file: file,
        };
        demod_struct
    }
}

impl<I> Iterator for FFTDemod<I> where I: Iterator<Item=Complex32> {
    type Item = f32;
    fn next(&mut self) -> Option<f32> {
        while self.window.len() < self.window_size {
            match self.iterator.next() {
                Some(value) => {
                    self.window.push_back(value);
                },
                None => return None
            }
        }
        for i in 0..self.window.len() {
            let pos: f32 = 0.5 + (i as f32) - (self.window_size as f32) / 2.;
            let gaussian_window: f32 = 3.;
            let window_function: f32 = (-pos * pos / (gaussian_window * gaussian_window)).exp();
            self.fft_input[i] = self.window[i] * window_function;
        }
        self.fft.process(& self.fft_input, &mut self.fft_output);
        self.window.pop_front();
        let mut count: f32 = 0.;
        let mut weighted: f32 = 0.;for n in 0..self.window_size  {
            let power = self.fft_output[n].norm();
            count += power;
            weighted += power * (n as f32 + 0.5);
        }
        match self.file {
            Some(ref mut file) => {
                for n in 0..self.window_size  {
                    let power = self.fft_output[n].norm();
                    count += power;
                    weighted += power * (n as f32 + 0.5);
                    //println!("{}", power);
                    write!(file, "{:.2},", self.fft_output[n].norm()).unwrap();
                }
                let treshold = self.window_size as f32 / 2.;
                let high = (treshold * 1.5) as i32;
                let low = (treshold * 0.5) as i32;
                write!(file, "{}\n", if weighted / count > treshold {high} else {low}).unwrap();
            },
            None => (),
        }
        Some(2. * (weighted / count) / (self.window_size as f32) - 1.)
    }
}

struct SymbolSync<T, I> {
    iterator: I,
    prev_sample: T,
    sample_symbol_ratio: f64,
    time_since_symbol: f64,
    sample_count: usize,
    file: Option<Box<Write>>,
}

impl<T, I> SymbolSync<T, I> where T: Eq + Copy + Default, I: Iterator<Item=T> {
    fn new(samples: I, ratio: f64, file: Option<Box<Write>>) -> SymbolSync<T, I> {
        let ss = SymbolSync {
            iterator: samples,
            prev_sample: T::default(),
            sample_symbol_ratio: ratio,
            time_since_symbol: 0.5,
            sample_count: 0,
            file: file,
        };
        ss
    }
}

impl<T, I> Iterator for SymbolSync<T, I> where T: Eq + Copy + Default + Display, I: Iterator<Item=T> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        loop {
            let new_sample = match self.iterator.next() {
                Some(value) => value,
                None => {
                    return None;
                }
            };
            match self.file {
                Some(ref mut file) => {
                    write!(file, "{},0\n", self.sample_count).unwrap();
                },
                None => (),
            }
            self.sample_count += 1;
            if self.prev_sample != new_sample {
                self.time_since_symbol = 0.5;
                self.prev_sample = new_sample;
            }
            self.time_since_symbol += self.sample_symbol_ratio;
            if self.time_since_symbol >= 1. {
                self.time_since_symbol -= 1.;
                match self.file {
                    Some(ref mut file) => {
                        write!(file, "{},{}\n", self.sample_count, new_sample).unwrap();
                        write!(file, "{},0\n", self.sample_count).unwrap();
                    },
                    None => (),
                }
                return Some(new_sample);
            }
        }
    }
}

struct IclickerDecode<T, I> where T: Copy + Eq + Zero, I: Iterator<Item=T> {
    iterator: I,
    window: u16,
}

impl<T, I> IclickerDecode<T, I> where T: Copy + Eq + Zero, I: Iterator<Item=T> {
    fn new(bits: I) -> IclickerDecode<T, I> {
        IclickerDecode {
            iterator: bits,
            window: 0,
        }
    }
}



impl<T, I> Iterator for IclickerDecode<T, I> where T: Copy + Eq + Zero + PartialOrd + std::fmt::Display, I: Iterator<Item=T> {
    type Item = IclickerPacket;
    fn next(&mut self) -> Option<IclickerPacket> {
        self.window = 0;
        while (self.window & 0xffff) != 0x5858 {
            let new_bit = match self.iterator.next() {
                Some(value) => {
                    value > T::zero()
                },
                None => {
                    return None;
                }
            };
            //println!("{:04x}", self.window);
            self.window = self.window << 1;
            if new_bit {
                self.window |= 0b1;
            }
        }
        let mut output: IclickerPacket = IclickerPacket::new();
        for i in 0..52 {
            let new_bit = match self.iterator.next() {
                Some(value) => {
                    value > T::zero()
                },
                None => {
                    return None;
                }
            };
            output.push_bit(new_bit);
        }
        Some(output)
    }
}

struct IclickerPacket {
    data: u64,
    len: usize,
}

impl IclickerPacket {
    fn new() -> IclickerPacket {
        IclickerPacket {
            data: 0,
            len: 0
        }
    }
    fn push_bit(&mut self, bit: bool) {
        if self.len > 64 {
            panic!("Cannot push more than 64 bits into packet.");
        }
        if bit {
            self.data |= 1 << (63 - self.len);
        }
        self.len += 1;
    }
    fn get_bit(&self, i: usize) -> u8 {
        return ((self.data >> (63 - i)) & 0x1) as u8;
    }
    fn get_range(&self, start: usize, len: usize) -> u64 {
        return (self.data << start) >> (64 - len);
    }
    fn bit_string(&self) -> String {
        let mut s = String::new();
        for i in 0..self.len {
            s.push_str(&self.get_bit(i).to_string());
            if i % 4 == 3 {
                s.push_str(" ");
            }
        }
        s
    }

    fn vote(&self) -> char {
        let vote_bits = self.get_range(40, 4);
        match vote_bits {
            0b0001 => 'A',
            0b0101 => 'B',
            0b1101 => 'C',
            0b1110 => 'D',
            0b1010 => 'E',
                 _ => '?'
        }
    }

    fn id(&self) -> u32 {
        let o: usize = 12;
        let b0: u8 = (
            (self.get_range(5 + o, 1) << 7) |
            (self.get_range(15 + o, 1) << 6) |
            (self.get_range(23 + o, 1) << 5) |
            (self.get_range(0 + o, 5) << 0)
        ) as u8;
        let b1: u8 = (
            (self.get_range(7 + o, 7) << 1) |
            (self.get_range(16 + o, 1))
        ) as u8;
        let b2: u8 = (
            (self.get_range(17 + o, 5) << 2) |
            (self.get_range(27 + o, 1) << 0)
        ) as u8;
        let b3: u8 = b0 ^ b1 ^ b2;
        return  ((b0 as u32) << 24) |
                ((b1 as u32) << 16) |
                ((b2 as u32) << 8) |
                ((b3 as u32) << 0);
    }
}

fn main() {
    let hw_frequency = 917e6;
    let sample_rate = 2.048e6;
    let bit_rate = 152.8e3;
    let process =
        Command::new("rtl_sdr")
        .arg("-f")
        .arg(hw_frequency.to_string())
        .arg("-s")
        .arg(sample_rate.to_string())
        .arg("-g")
        .arg("5.")
        .arg("-")
        .stdout(Stdio::piped())
        .stdin(Stdio::null())
        .spawn()
        .expect("failed to execute process");

    let mut demod_file = Box::new(BufWriter::new(File::create(Path::new("demod.csv")).unwrap()));
    let mut bits_file = Box::new(BufWriter::new(File::create(Path::new("bits.csv")).unwrap()));

    //let mut fft = rustfft::FFT::new(fft_size, false);
    let input = BufReader::new(process.stdout.unwrap());
    let mut stream = input.bytes()
        .map(|x| (x.unwrap() as f32 - 127.5) / 127.5)
        .tuples::<(_,_)>()
        .map(|x| Complex32::new(x.0, x.1));
    let mut stream = SquelchedSamples { iterator: stream, count: 0 };
    //let mut samples = Demod::new(samples, 32, Some(waterfall_file));
    let mut stream = QuadratureDemod::new(stream, Some(demod_file));
    let mut stream = SymbolSync::new(stream.map(|x| if x > 0. {1 as i8} else {-1 as i8}), bit_rate / sample_rate, Some(bits_file));
    let mut stream = IclickerDecode::new(stream);
    for y in stream {
        println!("vote: {}, id: {:08x}, packet contents: {}", y.vote(), y.id(), y.bit_string());
    }
}
