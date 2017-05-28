#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(unused_mut)]
#![allow(dead_code)]

extern crate core;
extern crate num;
extern crate itertools;
extern crate clap;

use std::io::{Write, Read, BufReader, BufWriter, self};
use std::process::{Command, Stdio};
use std::option;
use std::vec::Vec;
use std::boxed::Box;
use std::fmt::Display;
use std::path::Path;
use std::fs::File;
use std::str::FromStr;

use itertools::Itertools;

use num::complex::Complex32;
use num::complex::Complex;
use num::traits::Float;
use num::traits::Zero;
use core::ops::Neg;

use clap::{Arg, App, ArgGroup};

struct Treshold<I> {
    count: u8,
    iterator: I,
}

impl<I> Treshold<I> where I: Iterator<Item=Complex32> {
    fn treshold(&self, x: Complex32) -> bool {
        x.norm() > (0.1 as f32)
    }
}

impl<I> Iterator for Treshold<I> where I: Iterator<Item=Complex32> {
    type Item = Complex32;
    fn next(&mut self) -> Option<Complex32> {
        loop {
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
    window: u32,
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
        //0x858585 seems to be the sync word used by iclicker remotes
        while (self.window & 0xffffff) != 0x858585 {
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
        for i in 0..40 {
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
            if i % 8 == 7 {
                s.push_str(" ");
            }
        }
        s
    }

    fn vote(&self) -> char {
        let vote_bits = self.get_range(28, 4);
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
        //scrambled bytes
        let s0 = self.get_range(0,  8) as u8;
        let s1 = self.get_range(8,  8) as u8;
        let s2 = self.get_range(16, 8) as u8;
        let s3 = self.get_range(24, 8) as u8;
        //descrambled bytes
        let d0 = (s0 >> 3) | ((s2 & 1) << 5) | ((s1 & 1) << 6) | ((s0 & 4) << 5);
        let d1 = (s1 >> 1) | ((s0 & 1) << 7) | (s2 >> 7);
        //.......................................^.this one is questionable
        let d2 = ((s2 & 0x7c) << 1) | (s3 >> 5);
        //.........^..and these ones are not what Gourlay's paper shows
        let d3 = d0 ^ d1 ^ d2;
        //return as u32
        ((d0 as u32) << 24) | ((d1 as u32) << 16) | ((d2 as u32) << 8) | (d3 as u32)
    }

}

fn channel_to_frequency(value: &str) -> Option<f64> {
    let frequency = match value.to_uppercase().as_ref() {
        "AA" => 917, "AB" => 913, "AC" => 914, "AD" => 915,
        "BA" => 916, "BB" => 919, "BC" => 920, "BD" => 921,
        "CA" => 922, "CB" => 923, "CC" => 907, "CD" => 908,
        "DA" => 905, "DB" => 909, "DC" => 911, "DD" => 910,
        _ => 0,
    };
    if frequency > 0 { Some(frequency as f64 * 1e6) }
    else {
        match f64::from_str(value) {
            Ok(number) => Some(number),
            Err(e) => None
        }
    }
}

fn validate_channel(channel: String) -> Result<(), String> {
    match channel_to_frequency(channel.as_str()) {
        Some(frequency) => Ok(()),
        None => Err(format!("{} is not a valid channel.", channel))
    }
}

fn validate_gain(gain: String) -> Result<(), String> {
    match f64::from_str(gain.as_str()) {
        Ok(value) => {
            if value >= 0. {
                Ok(())
            }
            else {
                Err(format!("Gain must be nonnegative."))
            }
        }
        Err(e) => Err(format!("Gain must be a number."))
    }
}

fn validate_sample_rate(sample_rate: String) -> Result<(), String> {
    match f64::from_str(sample_rate.as_str()) {
        Ok(value) => {
            if value > 0. {
                Ok(())
            }
            else {
                Err(format!("Sample rate must be positive."))
            }
        }
        Err(e) => Err(format!("Sample rate must be a number."))
    }
}

fn validate_device_index(device_index: String) -> Result<(), String> {
    match usize::from_str(device_index.as_str()) {
        Ok(value) => {
            Ok(())
        }
        Err(e) => Err(format!("Device index must be a positive integer"))
    }
}

fn main() {
    let matches = App::new("Iclicker Sniffer")
        .version("0.1.0")
        .author("fb39ca4")
        .about("Receive iclicker packets with your RTL-SDR")
        .arg(Arg::with_name("channel")
            .short("c")
            .long("channel")
            .help("Sets the channel to use. Defaults to AA.")
            /*.possible_value("AA") .possible_value("AB") .possible_value("AC").possible_value("AD")
            .possible_value("BA") .possible_value("BB") .possible_value("BC").possible_value("BD")
            .possible_value("CA") .possible_value("CB") .possible_value("CC").possible_value("CD")
            .possible_value("DA") .possible_value("DB") .possible_value("DC").possible_value("DD")*/
            .validator(validate_channel)
            .takes_value(true))
        .arg(Arg::with_name("gain")
            .short("g")
            .long("gain")
            .help("RTL-SDR gain. Defaults to hardware automatic gain control.")
            .validator(validate_gain)
            .takes_value(true))
        .arg(Arg::with_name("sample_rate")
            .short("r")
            .long("sample_rate")
            .help("RTL-SDR sample rate. Defaults to 2.048 MHz.")
            .validator(validate_sample_rate)
            .takes_value(true))
        .arg(Arg::with_name("device")
            .short("d")
            .long("device")
            .help("RTL-SDR device index.")
            .validator(validate_device_index)
            .takes_value(true))
        .arg(Arg::with_name("stdin")
            .short("-s")
            .long("stdin")
            .help("Read data from stdin rather than from rtl_sdr")
            .conflicts_with_all(&["channel", "gain", "device_index"])
            .requires("sample_rate"))
        .arg(Arg::with_name("output")
            .short("-o")
            .long("output")
            .help("What to output")
            .possible_value("packets")
            .possible_value("bits")
            .default_value("packets")
            .takes_value(true))
        .get_matches();

    let input_stream : Box<Read>;
    let sample_rate;// = 2.048e6;
    if matches.is_present("stdin") {
        sample_rate = f64::from_str(matches.value_of("sample_rate").unwrap()).unwrap();
        input_stream = Box::new(io::stdin());
    }
    else {
        let frequency = channel_to_frequency(matches.value_of("channel").unwrap_or("AA")).unwrap();
        sample_rate = f64::from_str(matches.value_of("sample_rate").unwrap_or("900001")).unwrap();
        let mut command = Command::new("rtl_sdr");
        command.arg("-f").arg(frequency.to_string());
        command.arg("-s").arg(sample_rate.to_string());
        if let Some(gain) = matches.value_of("gain") {
            command.arg("-g").arg(gain);
        }

        let process = command
            .arg("-")
            .stdout(Stdio::piped())
            .stdin(Stdio::null())
            .spawn()
            .expect("failed to execute rtl_sdr process");

        input_stream = Box::new(process.stdout.unwrap());
    }

    let bit_rate = 152.8e3;

    let mut demod_file = Box::new(BufWriter::new(File::create(Path::new("demod.csv")).unwrap()));
    //let mut bits_file = Box::new(BufWriter::new(File::create(Path::new("bits.csv")).unwrap()));

    //let input = BufReader::new(process.stdout.unwrap());
    let mut stream = input_stream.bytes()
        .map(|x| (x.unwrap() as f32 - 127.5) / 127.5)
        .tuples::<(_,_)>()
        .map(|x| Complex32::new(x.0, x.1));
    let mut stream = Treshold { iterator: stream, count: 0 };
    let mut stream = QuadratureDemod::new(stream, Some(demod_file));
    let mut stream = SymbolSync::new(stream.map(|x| if x > 0. {1 as i8} else {-1 as i8}), bit_rate / sample_rate, None);

    if matches.value_of("output").unwrap() == "packets" {
        let mut stream = IclickerDecode::new(stream);
        for y in stream {
            println!("vote: {}, id: {:08x}, packet contents: {}", y.vote(), y.id(), y.bit_string());
        }
    }
    else {
        for y in stream {
            print!("{}", if y > 0 {1 as i8} else {0 as i8});
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            handle.flush().unwrap();
        }
    }
}
