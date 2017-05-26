# SmoothLifeJulia
SmoothLife implementation in Julia

Julia is a high level programming language which uses Phython and C as base. It can also use other C-based languages like C++ if necessary.
But unlike Python, Ruby, Matlab and other high level languages, Julia performs almost at C speed!
It also allows type-based functional programming, which helps to optimize code manually.

There seems to be no working implementation in Julia thus far, which might have to do with the a priori missing graphic support.
However, one can easily use and update Python images frames to do the output.

The idea has first been proposed in this paper: [Generalization of Conway's "Game of Life" to a continuous domain - SmoothLife](https://arxiv.org/abs/1111.1567)

## Note:
If you use Jupyter, the graphical output WILL NOT work. Jupyter will only output the last valid frame (at least in default configuration).

## TODO:
- [ ] stable universe
- [ ] working Jupyter output?
- [ ] code clean-up
- [ ] test cases
- [x] working graphics output
- [x] core