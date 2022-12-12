### A Pluto.jl notebook ###
# v0.19.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ fd79cce0-510f-11ec-2d74-fdc670099c49
begin
	using Pkg
	Pkg.activate(".")
	# Pkg.add("PlutoUI")
	using PlutoUI
	using Plots
	using Dates
	
	# global functions
	deg2rad(x) = x/180*pi
	rad2deg(x) = x/pi*180

	function get_md(J)
		date_origin = DateTime(2010, 1, 1, 0, 0, 0) + Day(J - 1)
		Dates.format(date_origin, "mm-dd")	
	end

	function ssh2time(ssh)
		w_hour = ssh/2
		w_sec = floor(Int, w_hour * 3600)
		time_begin = DateTime(2000, 1, 1, 12, 0, 0) - Dates.Second(w_sec)
		time_end = DateTime(2000, 1, 1, 12, 0, 0) + Dates.Second(w_sec)
		# time_begin, time_end
		time_begin, time_end = Dates.format(time_begin, "HH:MM:SS"), 
			Dates.format(time_end, "HH:MM:SS")
		time_begin, time_end	
	end
	
	# Solar Declination Angle （黄赤交角）
	function get_σ(J; to_deg = false)
		σ = 0.409 .* sin.(2*pi/365*J .- 1.39) # in [rad]
		if to_deg; σ = rad2deg.(σ); end
		σ
	end
end

# ╔═╡ 68cea5f1-4c4e-40af-8370-d2b7b616ab91


# ╔═╡ 0d8c5152-98f6-4925-a140-be892dadcd33
md"""
- 黄赤交角: σ
- 纬度: ϕ
- 时角: ω
"""

# ╔═╡ cc7a48d5-996b-435e-b2ca-8aa0ee847686
# # t = -12 : 12 # 夜间0 (-12)，第二天夜间为24(12)
# t = collect(0:24)
# t2 = t .- 12
# # t2 .* 15
# # ω = deg2rad(t2 .* 15)
# # length()
# # ω = deg2rad(t2 .* 15)
# # sinh = sin(ϕ) * sin(σ) .+ cos(ϕ) * cos(σ) * cos.(ω)
# # plot(sinh)
# # hline!([0], color = "black")

# ╔═╡ 148af5b3-aeba-400f-8d61-cc636732c631
begin
	# x = [1, 2, 3, 4]
	# clamp!(x, 2, 3)
	# x
end
# Longitude: 114.30898368551733
# Latitude: 30.585262910967806

# ╔═╡ b3f39acb-910d-463b-b174-c485378ebf2f
@bind J Slider(1:365)

# ╔═╡ c6f66283-0c2e-419d-a8c7-74b2dc428092
begin
	# d_r: reverse relative Earth-Sun distance 
	# I₀/ρ₀²
	function relativeDist(J)
		1 .+ 0.033 .* cos.(2*pi / 365 * J)
	end
	plot(relativeDist(1:365), title = "doy = $J", label = "", 
		framestyle = :box,
		ylab = "I₀/ρ₀²", xlab = "Doy of Year (DOY)")
	scatter!([J], [relativeDist(J)], markersize = 6, label = "")
end
# plot(collect(1:J))
# plot(x)

# ╔═╡ 72409fcc-45bd-4031-9c0f-9f2f16ad8e86
begin
	plot(get_σ(1:365; to_deg = true), 
		title = "doy = $J", label = "", 
		framestyle = :box,
		ylab = "Solar Declination Angle")
	scatter!([J], [get_σ(J; to_deg = true)], markersize = 6, label = "")
	# p
	# angle_SolarDecimation(1:10)
end

# ╔═╡ 32eae65b-cb25-421a-abc7-772fa28787b4
get_σ(J) |> rad2deg

# ╔═╡ 46f3563c-c9e1-4a33-8c8f-8bad7dad045a
J

# ╔═╡ 726c9e3b-f4f8-4534-a406-6ba28969a782
@bind ϕ_deg Slider(-90:90)

# ╔═╡ 2423c5b0-0dc7-42b1-8ead-33edcff0c2c6
begin
	# 太阳高度角sinh
	ω = collect(-2*pi:0.02:2*pi)
	ϕ = deg2rad(ϕ_deg)
	ϕ_deg
	# ϕ
	# # 获取日出时间和日落时间，即sinh = 0
	function get_ssh(ϕ, J)
		σ = get_σ.(J)
		# rad2deg(σ)
		tmp = - tan.(ϕ) .* tan.(σ)
		clamp!(tmp, -1, 1)
		# constrain in the range of [-1, 1]
		w = acos.(tmp)
		w_hour = w / pi * 12 # 距离中午的时间
		ssh = w_hour * 2
	end
	
	doy = 1:365
	time_sunrise = ssh2time(get_ssh(ϕ, [J])[1])
	title = "Sunshine Duration:
lat: $ϕ_deg, doy: $J, Year-Month: $(get_md(J)), 
sunrise and sunset = $time_sunrise"
	plot(doy, get_ssh(ϕ, doy), title = title, label = "")
	scatter!([J], get_ssh(ϕ, [J]), label = "")
	# get_ssh(ϕ, J)
end

# ╔═╡ Cell order:
# ╠═fd79cce0-510f-11ec-2d74-fdc670099c49
# ╠═68cea5f1-4c4e-40af-8370-d2b7b616ab91
# ╠═0d8c5152-98f6-4925-a140-be892dadcd33
# ╠═c6f66283-0c2e-419d-a8c7-74b2dc428092
# ╠═72409fcc-45bd-4031-9c0f-9f2f16ad8e86
# ╟─cc7a48d5-996b-435e-b2ca-8aa0ee847686
# ╟─2423c5b0-0dc7-42b1-8ead-33edcff0c2c6
# ╠═32eae65b-cb25-421a-abc7-772fa28787b4
# ╠═148af5b3-aeba-400f-8d61-cc636732c631
# ╠═46f3563c-c9e1-4a33-8c8f-8bad7dad045a
# ╠═b3f39acb-910d-463b-b174-c485378ebf2f
# ╠═726c9e3b-f4f8-4534-a406-6ba28969a782
