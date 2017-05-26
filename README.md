# SmoothLifeJulia
_SmoothLife implementation in Julia_

Julia is a high level programming language which uses Phython and C as base. It can also use other C-based languages like C++ if necessary.
But unlike Python, Ruby, Matlab and other high level languages, Julia performs almost at C speed!
It also allows type-based functional programming, which helps to optimize code manually.

There seems to be no working implementation in Julia _thus far_, which might have to do with the a priori _missing real-time rendering support_.
However, one can easily use and update **pyplot images** to make **real-time rendering**.

The idea has first been proposed in this paper: [Generalization of Conway's "Game of Life" to a continuous domain - SmoothLife](https://arxiv.org/abs/1111.1567)

## Note:
If you use Jupyter, real-time rendering **will not work**. Only the last frame will be shown in a frame. This is, unfortunately, an issues related to Jupyter (which I cannot bypass). 
_If someone has a workaround, let me know!_

## TODO:
- [ ] stable universe
- [ ] test cases
- [x] working graphics output
- [x] core
- [ ] code clean-up
- [ ] working Jupyter output?
