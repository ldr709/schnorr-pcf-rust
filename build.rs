fn main() {
    cc::Build::new()
        .cpp(true)
        .file("src/ntl_wrapper.cpp")
        .compile("ntl_wrapper");
    println!("cargo:rustc-link-lib=ntl");
}
