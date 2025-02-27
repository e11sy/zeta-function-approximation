run:
	julia --project=. src/main.jl

setup:
	julia --project=. -e 'import Pkg; Pkg.instantiate()'

test:
	julia --project=. test/runtests.jl
