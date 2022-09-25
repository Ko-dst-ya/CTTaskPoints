# Fixed Points.jl
module Points

import LinearAlgebra

export dist
export center
export neighbors
export Point
export Circle
export Square

# Следющие имена должны быть публичными:
# Point, neighbors, Circle, Square, center

"""
    Point(x, y)

Точка на декартовой плоскости.
"""
struct Point{T1<:Real, T2<:Real}
    x::T1
    y::T2
end

dist(p::Point) = sqrt(p.x^2 + p.y^2)

LinearAlgebra.norm(p::Point) = dist(p)
LinearAlgebra.dot(p1::Point, p2::Point) = p1.x * p2.x + p1.y * p2.y

Base.:+(p1::Point, p2::Point) = Point(p1.x + p2.x, p1.y + p2.y)
Base.:-(p1::Point, p2::Point) = Point(p1.x - p2.x, p1.y - p2.y)
Base.:-(p::Point) = Point(-p.x,-p.y)
Base.:*(α::Number, p::Point) = Point(α * p.x, α * p.y)
Base.:*(p::Point, α::Number) = α * p
Base.:/(p::Point, α::Number) = p * (1/α)

"""
    center(points) -> Point

Центр "масс" точек.
"""
function center(points)
	N = length(points)
	a = points[1]
	for i in 2:size(points, 1)
	     a += points[i]
	end
	return a/N
end

"""
    neighbors(points, origin, k) -> Vector{Point}

Поиск ближайших `k` соседей точки `origin` среди точек `points`.
"""
function neighbors(points, origin, k)
	k ≤ 0 && return []

	rho = [dist(p-origin) for p in points]

	D = Dict{Point, Real}()
	for i in 1:length(points)
    	D[points[i]] = rho[i]
    end
	ord = sort(collect(D), by=x->x[2])

	contcount=0
	res = []
	for i in 1:length(ord)
    	if isapprox(ord[i][2], 0.0)
        	contcount+=1
        	continue
		end
		if i<=k+contcount
			push!(res,ord[i][1])
		end
	end
	return(res)
end

"""
    Circle(o::Point, radius)

Круг с центром `o` и радиусом `radius`.
"""
struct Circle{TO<:Point, TR<:Real}
	o::TO
	radius::TR
end

Base.:in(p::Point, c::Circle) = dist(p-c.o) <= c.radius

"""
    Square(o::Point, side)

Квадрат с центром в `o` и стороной `side`. Стороны квадрата параллельны осям координат.
"""
struct Square{TO<:Point, TS<:Real}
	o::TO
	side::TS
end

Base.:in(p::Point, s::Square) = abs((s.o).x-p.x) <= s.side/2 && abs((s.o).y-p.y) <= s.side/2

"""
    center(points, area) -> Point

Центр масс точек `points`, принадлежащих области `area` (`Circle` или `Square`).	
"""
function center(points, area)
	return center(filter(p->p in area, points))
end

end # module
