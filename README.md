<h1 align="center">Ternimal</h1>
<h3 align="center">Simulate a lifeform in the terminal</h3>
<br>

Ter**n**i**m**al (note the spelling) is a program that draws an animated lifeform in the terminal using [Unicode block symbols](https://en.wikipedia.org/wiki/Block_Elements). It works in [most terminal emulators](#faq) and with most monospaced fonts.

![](https://user-images.githubusercontent.com/2702526/32404757-c4ee3230-c14e-11e7-9b5d-b48bd0fd2dab.gif)

From a practical perspective, the program is not very useful. It does, however, contain quite a bit of cool technology and math:

* "Glow" renderer capable of 1000+ frames per second (if the terminal can handle it)
* Dynamic generation of an everywhere differentiable movement path composed of circular arcs
* Skeletal deformation along the path
* [Fourier](https://en.wikipedia.org/wiki/Fourier_series)-based shape description allowing for many body forms to be realized and animated

Ternimal is also an exercise in minimalism and restraint. Written in just 1000 lines of Rust, it has *no dependencies* and consumes very few resources: 400 kB on disk, 3 MB of RAM and 4 % of a single CPU core with the default parameters. It implements its own linear algebra operations from scratch, as well as basic command line parsing and a simple random number generator.


## Building

Ternimal has no dependencies apart from the Rust Standard Library, and does not require Cargo for building. Only `rustc` (>= 1.20) must be installed, at which point Ternimal can be built with:

```
git clone https://github.com/p-e-w/ternimal.git
cd ternimal
rustc -O ternimal.rs
```


## Usage

Fundamentally, Ternimal does nothing more than color the distance field from a segmented spine moving along a meandering path. There are [many parameters](ternimal.rs#L88-L173) controlling this process, however, nearly all of which can be manipulated through the command line.

This makes the system very flexible. The following are just a few examples of what is possible:

### "Anaconda"

```
./ternimal length=100 segments=50 thickness=1,4,1,0 radius=6,12 gradient=0:#666600,0.5:#00ff00,1:#003300
```

![](https://user-images.githubusercontent.com/2702526/32404762-e5643794-c14e-11e7-81b2-bfa37809b128.gif)

Sine waves can be used to generate quite organic-looking shapes. In this case, a single half-wave forms the body of a snake.

### "Swarm"

```
./ternimal length=200 segments=50 thickness=0,4,19,0
```

![](https://user-images.githubusercontent.com/2702526/32404773-0e0a154c-c14f-11e7-8344-64e1d0e22617.gif)

Ternimal only renders a single model. However, thickness variations can give the appearance of multiple disconnected entities moving in a coordinated fashion.

### "Rainbow"

```
./ternimal length=20 thickness=70,15,0,5 padding=10 radius=5 gradient=0.03:#ffff00,0.15:#0000ff,0.3:#ff0000,0.5:#00ff00
```

![](https://user-images.githubusercontent.com/2702526/32404777-339d841a-c14f-11e7-97ee-b5f7a5ea87e3.gif)

The thickness function includes a time parameter. This makes it possible to define shape animations. Arbitrarily many Fourier series terms can be specified, enabling very complex animations.

### "Black Hole"

```
./ternimal speed=10 length=100 segments=5 thickness=13 gradient=0.5:#000000,0.8:#ffffff,1:#000000
```

![](https://user-images.githubusercontent.com/2702526/32404779-59c94552-c14f-11e7-8eea-d1786ff29e31.gif)

Creatively combining low segment counts, wide distance fields, and appropriately chosen gradients can produce output that seemingly goes beyond what is possible with spine-based rendering alone.


## FAQ

### What platforms and terminals are supported?

Ternimal has been tested on Linux, macOS, and Windows.

On Linux, almost all terminal emulators render Ternimal flawlessly. On macOS, iTerm2 or Alacritty are recommended. On Windows, PowerShell and WSL appear to work well, with ConEmu also working but lacking 24-bit colors.

### Why am I seeing random/weird colors?

Most likely because your terminal does not support 24-bit RGB color escape sequences. This in turn probably means that you are using macOS' default Terminal.app, which is the only major terminal emulator still missing that feature.

You have two options:

* Switch to a terminal that supports true color escape sequences. A well-researched list of such terminals can be found [here](https://gist.github.com/XVilka/8346728).
* Run Ternimal with the argument `true_color=false` to fall back to a 256-color palette, which is supported by practically every terminal emulator (but doesn't look as nice).

### Why am I seeing grid lines?

Either because your font's block characters do not completely fill the character cell, or because your terminal has line spacing greater than zero.

Ternimal works best with fonts that have a character aspect ratio as close to 2:1 as possible. The font in the screencasts is the wonderful [Iosevka](https://github.com/be5invis/Iosevka).

### Why does it look strange in the Linux (kernel) console?

The Linux console does not support Unicode fonts. It does, however, appear to recognize the Unicode *encoding*, and attempts to translate certain Unicode code points to code points in its internal encoding, which includes the block symbols required by Ternimal.

There appears to be a bug in this conversion, though. The *UPPER HALF BLOCK* character is translated correctly, but the *LOWER HALF BLOCK* is not, resulting in a striped pattern.

### How does it look in [cool-retro-term](https://github.com/Swordfish90/cool-retro-term)?

Like this:

![](https://user-images.githubusercontent.com/2702526/32404792-886a820e-c14f-11e7-9994-f27a0d048e39.gif)


## Contributing

Contributors are always welcome. However, **please file an issue describing what you intend to add before opening a pull request,** *especially* for new features! I have a clear vision of what I want (and do not want) Ternimal to be, so discussing potential additions might help you avoid duplication and wasted work.

By contributing, you agree to release your changes under the same license as the rest of the project (see below).


## License

Copyright &copy; 2017 Philipp Emanuel Weidmann (<pew@worldwidemann.com>)

Released under the terms of the [GNU General Public License, version 3](https://gnu.org/licenses/gpl.html)
