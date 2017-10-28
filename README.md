<h1 align="center">Ternimal</h1>
<h3 align="center">Simulate a lifeform in the terminal</h3>
<br>

Ter**n**i**m**al (note the spelling) is a program that draws an animated lifeform in the terminal using [Unicode block symbols](https://en.wikipedia.org/wiki/Block_Elements). It works in most terminal emulators and with most monospaced fonts.

From a practical perspective, the program is not very useful. It does, however, contain quite a bit of cool technology and math:

* "Glow" renderer capable of 1000+ frames per second (if the terminal can handle it)
* Dynamic generation of an everywhere differentiable movement path composed of circular arcs
* Skeletal deformation along the path
* [Fourier](https://en.wikipedia.org/wiki/Fourier_series)-based shape description allowing for many body forms to be realized and animated

Ternimal is also an exercise in minimalism and restraint. Written in just 1000 lines of Rust, it has *no dependencies* and consumes very few resources: 400 kB on disk, 3 MB of RAM and 3 % of a single CPU core with the default parameters. It implements its own linear algebra operations from scratch, as well as basic command line parsing and a simple random number generator.


## Building

Ternimal has no dependencies apart from the Rust Standard Library, and does not require Cargo for building. Only `rustc` (>= 1.20) must be installed, at which point Ternimal can be built with:

```
git clone https://github.com/p-e-w/ternimal.git
cd ternimal
rustc -O ternimal.rs
```


## Contributing

Contributors are always welcome. However, **please file an issue describing what you intend to add before opening a pull request,** *especially* for new features! I have a clear vision of what I want (and do not want) Ternimal to be, so discussing potential additions might help you avoid duplication and wasted work.

By contributing, you agree to release your changes under the same license as the rest of the project (see below).


## License

Copyright &copy; 2017 Philipp Emanuel Weidmann (<pew@worldwidemann.com>)

Released under the terms of the [GNU General Public License, version 3](https://gnu.org/licenses/gpl.html)
