# Not doing this the normal Julia way.

# Directly add path to the Julia LOAD_PATH if not already present
function add_to_Julia_path( p )
	if p ∉ LOAD_PATH
		push!( LOAD_PATH, p )
	end
end

add_to_Julia_path( "." )
add_to_Julia_path( "../git-submodules/julia-signals-systems/src/Sequences" )
add_to_Julia_path( "../git-submodules/julia-signals-systems/src/Processors" )
add_to_Julia_path( "../git-submodules/julia-signals-systems/src/ProcSeqs" )
add_to_Julia_path( "../git-submodules/julia-signals-systems/src/SNR" )
add_to_Julia_path( "../git-submodules/julia-signals-systems/src/DSPfns" )

add_to_Julia_path( "../git-submodules/julia-fixed-point" )

