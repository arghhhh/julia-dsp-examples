
# TODO: move this somewhere more appropriate

# DSP Components
# Major dependencies on both julia-fixed-point and julia-signals-systems

#--------------------------------------------------------------------------------

struct Truncate_lsbs_to{e} <: Processors.SampleProcessor
end
# value level:
function Processors.process( ::Truncate_lsbs_to{e1}, n, state=nothing ) where {e1}
        return FixPts.truncate_lsbs_to( e1, n ) , state
end
# type level:
function Base.eltype( ::Type{ Processors.Apply{I,Truncate_lsbs_to{e1}} }) where {I,e1} 
        return FixPts.truncate_lsbs_to( e1, eltype(I) )
end

#--------------------------------------------------------------------------------

struct Truncate_lsbs_by{e} <: Processors.SampleProcessor
end
# value level:
function Processors.process( ::Truncate_lsbs_by{e1}, n, state=nothing ) where {e1}
        return FixPts.truncate_lsbs_by( e1, n ) , state
end
# type level:
function Base.eltype( ::Type{ Processors.Apply{I,Truncate_lsbs_by{e1}} }) where {I,e1} 
        return FixPts.truncate_lsbs_by( e1, eltype(I) )
end

#--------------------------------------------------------------------------------

struct Round{T} <: Processors.SampleProcessor
end
# value level:
function Processors.process( ::Round{T}, n, state=nothing ) where {T}
        return Base.round( T, n ) , state
end
# type level:
function Base.eltype( ::Type{ Processors.Apply{I,Round{T}} }) where {I,T} 
        return T
end


#--------------------------------------------------------------------------------

# TODO: move image and alias somewhere more appropriate.

"Determine the frequency of the (i)th image of frequency f with sample rate fs"
function image( i, f, fs )
        @assert f <= fs/2 "input frequency should be in range 0..fs/2"
        add = iseven(i)
        mult = div( i+1, 2 )
        return mult * fs + (add ? f : -f )
end

"Determine what frequency f will alias to when sampled at fs"
function alias( f, fs )
        f1 = mod( f, fs )
        return f1 > fs/2 ? fs-f1 : f1
end
