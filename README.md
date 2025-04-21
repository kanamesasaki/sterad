# sterad

**sterad** is a collection of Rust crates for calculating radiation view factors and solid angles, both analytically and numerically. 
The main goal is to provide reliable, performant, and modular functionalities for scientific and engineering applications related to thermal radiation, nuclear engineering, illuminating engineering, computer graphics, and so on. 

## Overview

The **sterad** repository hosts multiple Rust crates, each specializing in different aspects of calculating solid angles, radiation view factors, and numerical (ray-tracing) solutions:

1. [sterad-view-factor](./crates/sterad-view-factor)  
   Analytical and semi-analytical computations of radiation view factors between surfaces.

2. [sterad-solid-angle](./crates/sterad-solid-angle)
   Analytical computations of solid angles for various geometries.

3. [sterad-ray-tracing](./crates/sterad-ray-tracing)  
   Calculating solid angles and view factors based on ray-tracing

## Features

- **Analytical Computations**: Compute solid angles and view factors for primitive shapes (rectangles, disks, spheres, cones, polygons, etc.).
- **Numerical Ray Tracing**: Monte Carlo-based ray tracing approaches for approximating view factors and solid angles.
- **Modular Design**: Each crate is independent and can be used separately or combined.
- **Rust Integration**: Written in pure Rust, designed to be efficient and safe.


## Installation

Each crate can be included in your `Cargo.toml`:

```toml
[dependencies]
sterad-view-factor = { git = "https://github.com/kanamesasaki/sterad", package = "sterad-view-factor" }
sterad-solid-angle = { git = "https://github.com/kanamesasaki/sterad", package = "sterad-solid-angle" }
sterad-ray-tracing = { git = "https://github.com/kanamesasaki/sterad", package = "sterad-ray-tracing" }
```

## Usage

Here are some examples of how to use the crates:

```rust
// Example for sterad-view-factor
use sterad_view_factor::diff_element_to_disk;

let h: f64 = 2.0;
let r: f64 = 1.0;
let vf: f64 = diff_element_to_disk::parallel_center(h, r).unwrap();
assert!((vf - 0.2).abs() < 1e-10);
```

```rust
// Example for sterad-solid-angle
use sterad_solid_angle::point_to_disk;
use std::f64::consts::PI;

let h = 4.0;
let r = 3.0;
let result = point_to_disk::center(h, r).unwrap();
assert!((result - 0.4 * PI).abs() < 1e-10);
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request if you have suggestions for improvements, bug fixes, or additional features.

## License

This project is licensed under the MIT License.
