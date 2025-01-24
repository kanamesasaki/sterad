# sterad

**sterad** is a collection of Rust crates for calculating radiation view factors and solid angles, both analytically and numerically. 
The main goal is to provide reliable, performant, and modular functionalities for scientific and engineering applications related to thermal radiation, nuclear engineering, illuminating engineering, computer graphics, and so on. 

## Overview

The **sterad** repository hosts multiple Rust crates, each specializing in different aspects of calculating solid angles, radiation view factors, and numerical (ray-tracing) solutions:

1. [sterad-numerical](./crates/sterad-numerical)  
   Ray-tracing and numerical methods for calculating solid angles, view factors, or other geometric measures.

2. [sterad-solid-angle](./crates/sterad-solid-angle)
   Analytical computations of solid angles for various geometries.

3. [sterad-view-factor](./crates/sterad-view-factor)  
   Analytical and semi-analytical computations of radiation view factors between surfaces.

## Features

- **Analytical Computations**: Compute solid angles and view factors for primitive shapes (rectangles, disks, spheres, cones, polygons, etc.).
- **Numerical Ray Tracing**: Monte Carlo-based ray tracing approaches for approximating view factors and solid angles.
- **Modular Design**: Each crate is independent and can be used separately or combined.
- **Rust Integration**: Written in pure Rust, designed to be efficient and safe.


## Installation

Each crate can be included in your `Cargo.toml`:

```toml
[dependencies]
sterad-numerical = { git = "https://github.com/kanamesasaki/sterad", package = "sterad-numerical" }
sterad-solid-angle = { git = "https://github.com/kanamesasaki/sterad", package = "sterad-solid-angle" }
sterad-view-factor = { git = "https://github.com/kanamesasaki/sterad", package = "sterad-view-factor" }
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request if you have suggestions for improvements, bug fixes, or additional features.

## License

This project is licensed under the MIT License.
