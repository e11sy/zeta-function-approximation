run:
	julia --project=. src/main.jl

vis:
	julia --project=. src/visualize_n.jl

setup:
	julia --project=. -e 'import Pkg; Pkg.instantiate()'

test:
	julia --project=. test/runtests.jl
