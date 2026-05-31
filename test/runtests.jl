using Test

# Runtime stubs for translated routines
numtasks = 1
pid = 0
status = nothing
ierr = 0
MPI_wait(args...) = nothing
sendborders(args...) = (0, 0, 0, 0)
sendbordersflood(args...) = (0, 0, 0, 0)
Ksat(::Integer) = 10.0
klatfactor(::Integer) = 1.0
const pi4 = 4π

include("../src/Groundwater/flowdir.jl")
include("../src/Groundwater/GW2river.jl")
include("../src/Groundwater/moveqrf.jl")
include("../src/Groundwater/lateralflow.jl")
include("../src/Groundwater/flooding.jl")
include("../src/Groundwater/rivers_kw_flood.jl")

@testset "flowdir translation" begin
  fd = zeros(Int, 3, 3)
  fd[2, 2] = 1
  @test flowdir(3, 1, 3, fd, 2, 2) == (3, 2)

  fd[2, 2] = 8
  @test flowdir(3, 1, 3, fd, 2, 2) == (1, 1)
end

@testset "gw2river branches" begin
  imax, js, je, nzg = 3, 1, 3, 1
  slz = [-1.0, -2.0]
  deltat = 1.0
  soiltxt = ones(Int, 2, imax, je)
  landmask = zeros(Int, imax, je)
  landmask[2, 2] = 1
  maxdepth = fill(2.0, imax, je)
  riverdepth = fill(1.0, imax, je)
  width = zeros(Float64, imax, je)
  width[2, 2] = 2.0
  length = fill(3.0, imax, je)
  area = fill(6.0, imax, je)
  fdepth = fill(2.0, imax, je)
  qrf = zeros(Float64, imax, je)

  wtd = fill(-0.5, imax, je)
  gw2river(imax, js, je, nzg, slz, deltat, soiltxt, landmask, wtd, maxdepth, riverdepth, width, length, area, fdepth, qrf)
  @test qrf[2, 2] > 0.0

  wtd[2, 2] = -1.5
  gw2river(imax, js, je, nzg, slz, deltat, soiltxt, landmask, wtd, maxdepth, riverdepth, width, length, area, fdepth, qrf)
  @test qrf[2, 2] ≈ -1.0 atol = 1e-8

  wtd[2, 2] = -2.5
  gw2river(imax, js, je, nzg, slz, deltat, soiltxt, landmask, wtd, maxdepth, riverdepth, width, length, area, fdepth, qrf)
  @test qrf[2, 2] ≈ -1.0 atol = 1e-8
end

@testset "moveqrf routing" begin
  fd = zeros(Int, 3, 3)
  fd[2, 2] = 1
  qrf = zeros(Float64, 3, 3)
  qrf[2, 2] = 2.0
  area = ones(Float64, 3, 3)
  width = ones(Float64, 3, 3)
  width[2, 2] = 0.5

  moveqrf(3, 1, 3, fd, qrf, area, width)
  @test qrf[2, 2] == 0.0
  @test qrf[3, 2] ≈ 2.0 atol = 1e-8
end

@testset "lateralflow4 flat head" begin
  wtd = fill(-1.0, 3, 3)
  qlat = zeros(Float64, 3, 3)
  fdepth = fill(1.0, 3, 3)
  topo = fill(100.0, 3, 3)
  landmask = ones(Int, 3, 3)
  area = ones(Float64, 3, 3)
  klat = ones(Float64, 3, 3)
  xlat = zeros(Float64, 3, 3)

  lateralflow4(3, 3, 1, 3, wtd, qlat, fdepth, topo, landmask, 1.0, area, klat, xlat, 1.0)
  @test qlat[2, 2] == 0.0
end

@testset "flooding simple redistribution" begin
  fd = zeros(Int, 3, 3)
  fd[2, 2] = 1
  bfd = zeros(Int, 3, 3)
  topo = zeros(Float64, 3, 3)
  area = ones(Float64, 3, 3)
  riverwidth = zeros(Float64, 3, 3)
  riverlength = ones(Float64, 3, 3)
  riverdepth = ones(Float64, 3, 3)
  floodheight = zeros(Float64, 3, 3)
  floodheight[2, 2] = 1.0
  delsfcwat = zeros(Float64, 3, 3)

  flooding(3, 1, 3, 1.0, fd, bfd, topo, area, riverwidth, riverlength, riverdepth, floodheight, delsfcwat)
  @test delsfcwat[2, 2] < 0.0
  @test sum(delsfcwat) ≈ 0.0 atol = 1e-8
end

@testset "rivers_kw_flood smoke test" begin
  imax, js, je = 3, 1, 3
  deltat, dtlr = 1.0, 1.0
  fd = zeros(Int, imax, je)
  bfd = zeros(Int, imax, je)
  fd[2, 2] = 1
  qnew = zeros(Float64, imax, je)
  qnew[2, 2] = 1.0
  qs = zeros(Float64, imax, je)
  qrf = zeros(Float64, imax, je)
  delsfcwat = zeros(Float64, imax, je)
  slope = fill(0.01, imax, je)
  depth = fill(1.0, imax, je)
  width = fill(2.0, imax, je)
  length = fill(10.0, imax, je)
  maxdepth = fill(1.0, imax, je)
  area = fill(10.0, imax, je)
  riverarea = fill(20.0, imax, je)
  floodarea = fill(0.0, imax, je)
  riverchannel = fill(20.0, imax, je)
  qmean = zeros(Float64, imax, je)
  floodheight = zeros(Float64, imax, je)
  topo = zeros(Float64, imax, je)

  rivers_kw_flood(imax, js, je, deltat, dtlr, fd, bfd, qnew, qs, qrf, delsfcwat,
                  slope, depth, width, length, maxdepth, area, riverarea, floodarea,
                  riverchannel, qmean, floodheight, topo)

  @test qnew[2, 2] > 0.0
  @test qmean[2, 2] > 0.0
end
