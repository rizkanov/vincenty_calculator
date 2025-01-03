# Rust Project vincenty-calculator
This project was created to compile a Rust project and display the resulting WebAssembly (WASM) on a web page. Both vincety_calculator and fe-vincenty-calculator need to be stored in the same directory to run in development mode.

Clone project

## Generate Rust to Web

```bash
wasm-pack build --target web
```

## Make package available to npm

```bash
wasm-pack build --target bundler
```
