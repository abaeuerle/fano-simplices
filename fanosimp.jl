###########################################################
###   A JULIA PACKAGE FOR CLASSIFYING FANO SIMPLICES    ###
###########################################################

using Combinatorics, Primes, AbstractAlgebra, Base.Threads

global count = Atomic{Int}(0)

############################
###   Helper functions   ###
############################

struct fanosimp
    pl::Vector{Vector{Int}}
    ul::Vector{Vector{Rational{Int}}}
end

function wstring( w::Vector{Int} )
    s = filter( x -> !isspace( x ), string( w ) )
    s = replace( s, "["=>"" )
    s = replace( s, "]"=>"" )
    s = replace( s, ","=>"-" )
    return s
end 

# fillvec!( v, n ) extends integer list v of length m by n-m zeros.
function fillvec!( v::Vector{Int}, n::Int )
    append!( v, [ zero(Int) for i = 1:n-length(v) ] )
end

# divisors( n ) returns list of all positive divisors of n.
# If n = 0, returns empty list.
function divisors( n::Int )
    if n < 0
        return divisors( -n )
    end
    if n == 0
        return Int[]
    end
    if isprime(n)
        return [one(Int), n]
    end 
    f = Primes.factor(n)
    d = [one(Int)]
    for (k, v) in f
        c = [k^i for i in 0:v]
        d = d*c'
        d = reshape(d, length(d))
    end
    return sort!(d)
end

function str2vec( s::String )
    return parse.( Int, split( s[2:end-1], "," ) )
end

function str2vecvec( s::String )
    sl = split( replace(s[3:end-2], " "=>"" ), "],[" )
    return [parse.(Int, split(s, ",")) for s in sl]
end

function swap!( l, i, j )
    t = l[i]
    l[i] = l[j]
    l[j] = t
end

function basechange( z::Int, b::Int, n::Int )
    zl = Int[]
    for j = 1:n
        t = sum( [ b^(k-1)*zl[k] for k = 1:length( zl ) ] )
        append!( zl, div( (z - t) % b^j, b^(j-1) ) )
    end
    return zl
end


######################################################
###   Unit fraction partitions and weight systems  ###
######################################################

function isintegral( q::Rational{Int} )
    return denominator( q ) == one(Int)
end

function isunitfraction( q::Rational{Int} )
    return isintegral( one(Int)//q )
end

function isufp( ufl::Vector{Rational{Int}} ) 
    return isunitfraction( sum( ufl ) )
end

# allufp( q, n ) computes all sorted unit fraction partitions of q of length n.
# Ie. all n-tuples [1//a[1],...,1//a[n]] with a[1] <= ... <= a[n],
# such that q = 1//a[1] + ... + 1//a[n].
function allufp( q::Rational{Int}, n::Int )
    if q <= 0 || n <= 0
        return []
    end
    if n == 1
        if isunitfraction( q )
            return [[q]]
        end
        return []
    end
    ufl = []
    for j = ceil( Integer, 1/q ) : floor( Integer, n/q )
        ufltemp = allufp( q - 1//j, n-1 )
        append!( ufl, sort( append!( ufd, 1//j), rev = true ) for ufd in ufltemp )
    end
    return union( ufl )
end

# ufp2w( A ) takes a unit fraction partition A = [1//a[1],...,1//a[n]] and returns
# the corresponding reduced weight system Q(A).
function ufp2w( ufl::Vector{Rational{Int}} )
    t = lcm( denominator.(ufl) )
    return numerator.(t.*ufl)
end

# w2ufp( gi, w ) takes an integral weight system w = [w[1],...,w[n]] and returns
# the cooresponding unit fraction partition A(w) of 1//g in the form [1//a[1],...,1//a[n]].
# This only works if g is a multiple of the index of w.
function w2ufp( gi::Int, w::Vector{Int} )
    ww = sum(w)
    return w.//(gi*ww)
end

# iswellformed( w ) checks whether the integral weight system w is well-formed.
function iswellformed( w::Vector{Int} )
    for i = 1 : length(w)
        if gcd( deleteat!( copy( w ), i) ) != 1
            return false
        end
    end
    return true
end


#####################################
###   Procedures for P-matrices   ###
#####################################

# We store integral matrices P as lists of lists of integers pl,
# where the inner lists are the COLUMNS of P. Ie. the matrix
#
#     | a11   ...   a1n |
# P = | .           .   |
#     | .           .   |
#     | am1   ...   amn |
#
# is stored as pl = [[a11,...,am1],...,[a1n,...,amn]].

# IN WHAT FOLLOWS, P IS ALWAY ASSUMED TO BE A d x (d+1) INTEGER MATRIX
# WHOSE COLUMNS ARE THE VERTICES OF A FANO SIMPLEX \delta. (We just call them P-matrices)
# pl IS THUS THE LIST OF VERTICES OF \delta, STORED AS LISTS OF INTEGERS.

# gorform( pl, w, i ) computes the i-th Gorenstein form of the simplex \Delta = \Delta(pl).
# - pl must be a list of d+1 list of d integers, being the vertices of a d-dimensional Fano simplex.
# - w is the weight system of \Delta.
# - 1 <= i <= d+1.
function gorform( pl::Vector{Vector{Int}}, w::Vector{Int}, i::Integer )::Vector{Rational{Int}}
    n = length(pl) - 1
    wZ = sum( w )
    ul = []
    for j = 1:n
        append!( ul, j == i ? qval( pl, ul ) + wZ//( w[i] * pl[i][i] ) : qval( pl, ul ) )
    end
    return ul
end

# qval( pl, ql ) inductively computes the rational numbers
# q[k] from the proof of Proposition 5.2.
# - pl is the list of vertices of a Fano simplex
# - ql = [q[1],...,q[k]-1]
function qval( pl::Vector{Vector{Int}}, ql )
    n = length(ql)+1
    if n == 1
        return -one(Int)//one(Int)
    end
    return -(1 + sum( pl[n][j] * ql[j] for j = 1:n-1 ) )//pl[n][n]
end

# uval( pl, w, i, ul ) returns the next entry of the i-th Gorenstein form
# of the d-dimensional Fano simplex \Delta = \Delta(pl).
# - pl is the list of vertices of \Delta.
# - w is the weight system of \Delta.
# - 1 <= i <= d+1.
# - ul is the list of the previous entries of the i-th Gorenstein form.
#   Ie. if ul = [u[i,1],...,u[i,k]] then uval returns u[i,k+1].
function uval( pl::Vector{Vector{Int}}, w::Vector{Int}, i::Int, ul )
    n = length(ul)+1
    wZ = sum( w )
    return n == i ? qval( pl, ul ) + wZ//( w[i] * pl[i][i] ) : qval( pl, ul )
end

# qvcol does the same as qval, where cl is the n-th column of pl.
function qvcol( cl::Vector{Int}, ql )
    n = length(ql)+1
    if n == 1
        return - 1//1
    end
    return - (1 + sum( cl[j] * ql[j] for j = 1:n-1 ) )//cl[n]
end

# qvcol does the same as uval, where cl is the n-th column of pl.
function uvcol( cl::Vector{Int}, w::Vector{Int}, i::Int, ul )
    n = length(ul)+1
    wZ = sum( w )
    return n == i ? qvcol( cl, ul ) + wZ//( w[i] * cl[i] ) : qvcol( cl, ul )
end

# isintegral( v ) checks whether the entries of the rational vector v are integers.
function isintegral( v::Vector{Rational{Int}} )
    return all( isintegral.(v) )
end

# isprimitive( v ) checks whether the rational vector v is a primitive lattice point.
function isprimitive( v::Vector{Rational{T}} ) where {T<:Integer}
    prod( [ denominator(q) for q in v ] ) == one(T) || return false
    return gcd( [ numerator(q) for q in v ] ) == one(T)
end

# pl2lgi( pl, w ) returns the local Gorenstein indices [g[1],...,g[d+1]] of the
# d-dimensional Fano simplex \Delta = \Delta(pl).
# - pl is the list of vertices of \Delta.
# - w is the weight system of \Delta.
function pl2lgi( pl::Vector{Vector{Int}}, w::Vector{Int} )
    lgi = Int[]
    for i = 1:length( pl )
        ui = gorform( pl, w, i )
        append!( lgi, lcm( [ denominator( q ) for q in ui ] ) )
    end
    return lgi
end

# hnfl(pl) computes the Hermite normal form H of the matrix P
# corresponding to pl and returns it as a list of the columns of H.
# - pl is a list of the columns of P.
function hnfl( pl::Vector{Vector{Int}} )
    n = length( pl )
    m = matrix( ZZ, reduce( hcat, pl ) )
    h = hnf( m )
    return [ [ convert( Int, h[i,j] ) for i = 1:n-1 ] for j = 1:n ]
end

# lexmin( pllist ) returns the lexicographic minimum of pllist.
# - pllist is a list of pl-list, ie. of lists of lists of integers.
function lexmin( pllist::Vector{Vector{Vector{Int}}} )
    m = pllist[1]
    n1 = length(m)
    n2 = length(m[1])
    for i = 2:length(pllist)
        l1 = [ m[i][j] for j = 1:n2 for i = 1:n1 ]
        pl = pllist[i]
        l2 = [ pl[i][j] for j = 1:n2 for i = 1:n1 ]
        if l2 < l1
            m = pl
        end
    end
    return m
end

# normalform!( pl, w ) replaces pl by its normal form as defined in Definition 5.3.
# - pl is a list of the vertices of a Fano simplex \Delta.
# - w is the weight system of \Delta.
function normalform!( pl::Vector{Vector{Int}}, w::Vector{Int} )
    n = length( pl )
    lgi = pl2lgi( pl, w )
    for i = 1:n-1
        for j = i+1:n
            if w[j] > w[i]
                swap!( w, i, j )
                swap!( lgi, i, j )
                swap!( pl, i, j )
            elseif w[i] == w[j]
                if lgi[j] > lgi[i]
                    swap!( lgi, i, j )
                    swap!( pl, i, j )
                end
            end
        end
    end
    indexgroup = [[1]]
    for i = 2:n
        if w[i-1] == w[i] && lgi[i-1] == lgi[i]
            append!(last(indexgroup),i)
        else
            push!( indexgroup, [i] )
        end
    end
    plpermute = [Vector{Int}[]]
    for g in indexgroup
        plltemp = Vector{Vector{Int}}[]
        gperms = permutations(g)
        for gp in gperms
            for pp in plpermute
                push!( plltemp, append!( copy(pp), [ pl[i] for i in gp ] ) )
            end
        end
        plpermute = copy( plltemp )
    end
    hnflist = [ hnfl( p ) for p in plpermute ]
    pl[:] = lexmin( hnflist )
    return pl
end

# weights( pl ) computes the weight system of the Fano simplex \Delta = \Delta(pl)
# - pl is the list of vertices of \Delta.
function weights( pl::Vector{Vector{Int}} )
    w = Rational{Int}[]
    for i in eachindex(pl)
        m = reduce( hcat, deleteat!( copy(pl), i ) )
        append!( w, Rational{Int}( abs( det( matrix( QQ, m ) ) ) ) )
    end
    return w
end

# multiplicity( pl ) computes the multiplicity of the Fano simplex \Delta = \Delta(pl)
# - pl is the list of vertices of \Delta.
function multiplicity( pl::Vector{Vector{Int}} )
    w = weights( pl )
    d = lcm( denominator.(w) )
    return gcd( numerator.(d.*w) )//d
end


##########################################
###   Procedures for classification    ###
##########################################

# dcol( plaffine ) computes the last column of the P-matrix of \Delta, assuming
# plaffine is the list of the first d columns.
# - w is the weight system of \Delta.
function dcol( plaffine::Vector{Vector{Int}}, w::Vector{Int} )
    n = length( w )-1
    wl = last(w)
    return [ -one(Int)//wl * sum( [ plaffine[i][j] * w[i] for i = 1:n ] ) for j = 1:n ]
end

# affinecandidates( dmaxl, w, lgi ) returns a list of all candidates for the first
# d columns of a P-matrix in hermite normal form of a Fano simplex \Delta with
# weight system w and local Gorenstein indices lgi.
# The i-th diagonal entry of P is a divisor of the (i-1)-th entry of dmaxl, i >= 2. 
function affinecandidates( dmaxl::Vector{Int}, w::Vector{Int}, lgi::Vector{Int} )
    outputlist = fanosimp[]
    n = length( w )
    m = length( dmaxl )
    dvals = divisors( last( dmaxl ) )
    if m == 1
        ul = [ [ uvcol( [1], w, i, [] ) ] for i=1:n ]
        for d in dvals
            for a = 0:d-1
                pl = [ [1], [a,d] ]
                ul2 = [ push!( copy( ul[i] ), uvcol( last( pl ), w, i, ul[i] ) ) for i=1:n ]
                if gcd( a, d ) == 1 && all( isintegral.( lgi .* ul2 ) )
                    push!( outputlist, fanosimp( pl, ul2 ) )
                end
            end
        end
        return outputlist
    end
    pllt = affinecandidates( deleteat!( copy(dmaxl), length(dmaxl) ), w, lgi )
    for l in pllt
        pl = l.pl
        ul = l.ul
        for d in dvals
            iterarray = reshape( 1:d^m, Tuple( d for i=1:m ) )
            cartiter = CartesianIndices(iterarray)
            for cartindex in cartiter
                offdiagonalvals = collect( Tuple(cartindex) )
                offdiagonalvals .-= 1
                newcolumn = append!( offdiagonalvals, d )
                if gcd( newcolumn ) == 1
                    ullast = [ uvcol( newcolumn, w, i, ul[i] ) for i=1:n ]
                    if all( isintegral( lgi .* ullast ) )
                        pl2 = push!( copy( pl ), newcolumn )
                        ul2 = [ push!( copy( ul[i] ), ullast[i] ) for i=1:n ]
                        push!( outputlist, fanosimp( pl2, ul2 ) )
                    end
                end
            end
        end
    end
    return outputlist
end

# classify_fanosimp( g, w ) classifies all fano simplices of Gorenstein index g
# with weight system w. g must be a multiple of the index of w.
function classify_fanosimp( gi::Int, w::Vector{Int} )
    outputlist = []
    wZ = sum( w )
    n = length( w ) - 1
    ufp_a = denominator.( w2ufp( gi, w ) )
    dmaxl = [ gcd(ufp_a[1],ufp_a[2]) ]
    lambda_max = numerator(prod(ufp_a)//lcm(ufp_a)^2)
    for k = 3:n
        push!( dmaxl,gcd( ufp_a[k], last(w)*lambda_max ) )
    end
    affcands = affinecandidates( dmaxl, w, [gi for i = 1:length(w)] )
    templist = []
    for plaff in affcands
        ul = plaff.ul
        lgi = [ lcm( denominator.(u) ) for u in ul ]
        if lcm( lgi ) == gi
            plcols = [ fillvec!( copy(v), n ) for v in plaff.pl ]
            dc = dcol( plcols, w )
            if isprimitive( dc )
                push!( plcols, numerator.(dc) )
                push!( templist, plcols )
            end
        end
    end
    for pl in templist
        normalform!( pl, w )
        if !( pl in outputlist )
            push!( outputlist, pl )
        end
    end
    return outputlist
end

# fullclass_fanosimp( d, g ) performs the full classification of d-dimensional Fano
# simplices of Gorenstein index g.
# Prints information on the current progress of the classification.
# The function creates a folder called "temp" in the current working directory.
# The classifaction results are written as text files for each weight vector in the temp-folder.
# The function is multithreaded.
function fullclass_fanosimp( n::Int, gi::Int )
    mkpath("./temp")
    println( "Classification of Fano simplices of dimension " * string( n ) * " and Gorenstein index " * string( gi ) * "." )
    println( "Using " * string( nthreads() ) * " threads." )
    println( "Computing list of weight systems..." )
    # compute all ufp's of g of length d+1.
    ufplist = allufp( 1//gi, n+1 )
    # compute weight systems and filter by well-formedness.
    wlist = ufp2w.( ufplist )
    wlist = filter( w -> iswellformed(w), wlist )
    println("Done. List of " * string( length( wlist ) ) * " well-formed weight systems.")
    println( "Starting classification. This might take a while." )
    # Perform classification:
    global count = Atomic{Int}(0)
    @threads for w in wlist
        wrapper( gi, w, length( wlist ) )
    end
    println( "Done." )
    # Count number of Fano simplices:
    plcount = 0
    for w in wlist
        inputfile = "./temp/fanosimp-d=" * string( n ) * "-gi=" * string( gi ) * "-w=" * wstring( w ) * ".txt"
        open( inputfile, "r" ) do fin
            lines = readlines( fin )
            plcount += length( lines )
        end
    end
    println( "Found " * string( plcount ) * " Fano simplices of dimension " * string(n) * " and Gorenstein index " * string(gi) * "." )
end

# Wrapper for the subclassification. Prints the current status.
# Writes subclassification to a text file. Needs folder "temp" in working directory.
function wrapper( gi::Int, w::Vector{Int}, nw::Int )
    n = length( w ) - 1
    tempfile = "./temp/fanosimp-d=" * string( n ) * "-gi=" * string( gi ) * "-w=" * wstring( w ) * ".txt"
    println( string( w ) * " computing..." )
    simplist = classify_fanosimp( gi, w )
    atomic_add!( count, 1 )
    println( string( w ) * " done. " * string( count[] ) * "/" * string( nw ) )
    open( tempfile, "w" ) do file
        for pl in simplist
            write( file, string( pl ) * "\n" )
        end
    end
end
