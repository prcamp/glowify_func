using PointsNCurves
#const LineSamplePoints = StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
struct DiscreteDomain
    Nx::Int
    Ny::Int
    minmin::Point
    maxmax::Point
end

function draw_glowified_fun(filename, horz_planks, vert_planks, unit)
  DrawSVGCurvesNoScaling(filename, map(c -> Curve(c), vcat(horz_planks,vert_planks)),unit)
end


function glowify_curves_on_grid(vert_curves::Curves, horz_curves::Curves, mat_width::Float64,x_length::Float64, y_length::Float64, z_length::Float64, z_min::Float64)
  # vert_curves will get cut from above
  # horz_curves will get cut from below

  # Assuming even spacing of grid lines.
  # I have in mind that the minimum and maximum x coordinates are uniform over vert_curves and horz_curves
  # Major assumption: the "x" coordinates of the points in vert_curves and horz_curves are monotonically non-decreasing!
  
  depth_factor = .5
  cut_width = mat_width


  Nvert = length(vert_curves)
  Nhorz = length(horz_curves)

  # we'll scale the values according to z_length and the input min, max of the function
  all_func_max = max(mapreduce(c -> maxy(c),max,vert_curves),mapreduce(c -> maxy(c), max, horz_curves))
  all_func_min = min(mapreduce(c -> miny(c),min,vert_curves),mapreduce(c -> miny(c), min, horz_curves))

  # scale the "y" coordinates of the curves by z_length, and shift by z_min:
  scale_z(c::Curve) = map!(p -> Point(p.x, z_min + z_length*((p.y - all_func_min)/(all_func_max - all_func_min))),c)

  for c in vert_curves
    scale_z(c)
  end
  for c in horz_curves
    scale_z(c)
  end
  
  # Now the depth axis is scaled as desired.

  # (re)compute the intersection points given the curves.
  horz_input_min = mapreduce(c -> minx(c), min, horz_curves)
  horz_input_max = mapreduce(c -> maxx(c), max, horz_curves)
  vert_input_min = mapreduce(c -> minx(c), min, vert_curves)
  vert_input_max = mapreduce(c -> maxx(c), max, vert_curves)

  horz_disc_vals = linspace(horz_input_min, horz_input_max,  Nvert + 2)[2:end-1]
  vert_disc_vals = linspace(vert_input_min, vert_input_max,  Nhorz + 2)[2:end-1]

  # Now compute the (minimum) of the depths of the horz_curves and vert_curves where they should intersect:
  horz_vals_along_grid = Array{Float64}(Nvert,Nhorz)
  vert_vals_along_grid = Array{Float64}(Nvert,Nhorz)
  cut_depth_from_bottom = Array{Float64}(Nvert,Nhorz)

  box_curve(c::Curve, xmin, xmax)::ClosedCurve = ClosedCurve(vcat([Point(xmin,0)],c, [Point(xmax,0)]))
  horz_planks = map(c -> box_curve(c,horz_input_min,horz_input_max),horz_curves)
  vert_planks = map(c -> box_curve(c,vert_input_min,vert_input_max),vert_curves)
  
  horz_cuts = ClosedCurves()
  vert_cuts = ClosedCurves()
  
    # from above
  vert_cutout(x_coord,cut_depth) = ClosedCurve(
    [Point(x_coord+cut_width/2,z_length+z_min+5),
    Point(x_coord+cut_width/2,cut_depth),
    Point(x_coord-cut_width/2,cut_depth),
    Point(x_coord-cut_width/2,z_length+z_min+5)            
    ])
    # from below
  horz_cutout(x_coord,cut_depth) = ClosedCurve(
    [
    Point(x_coord-cut_width/2,-5),
    Point(x_coord+cut_width/2,-5),                      
    Point(x_coord+cut_width/2,cut_depth),            
    Point(x_coord-cut_width/2,cut_depth)                        
    ])

  simple_cut(p1,p2,cut_depth) = Curve(
    [
    p1,
    Point(p1.x,cut_depth),
    Point(p2.x,cut_depth),
    p2
    ]
    )

  for (i,horz) in enumerate(horz_curves) # Nhorz of these. They cross the vertical curves
    #cur_horz_idx = 1
    for (j,vert) in enumerate(vert_curves) # Nvert of these. They cross the horizontal curves
    #  cur_vert_idx = 1
      horz_coord = horz_disc_vals[j]
      vert_coord = vert_disc_vals[i]
      # Do horz or verts have values at those coordinates?
      # increment the horizontal and vertial indices until the "x" coordinates 
#       while (horz[cur_horz_idx].x < horz_coord)
#         cur_horz_idx += 1
#       end
#       while (vert[cur_vert_idx].x < vert_coord)
#         cur_vert_idx += 1
#       end
#       min_horz_func_val = horz[cur_horz_idx].y
#       min_vert_func_val = vert[cur_vert_idx].y            
    horzfunccut = vert_cutout(horz_coord, z_min)
    ~,~,horz_func_vals_at_intersection = masked_occluded_intersections(horz_planks[i],horzfunccut)
    # if (isempty(horz_func_vals_at_intersection))
    #   println("length of horz_plank is $(length(vert_planks[j]))")
    #   #PointsNCurves.DrawSVGCurves("vert_curves_$j.svg",map(c -> Curve(c), [vert_planks[j]]))
    # else
    #   min_vert_func_val = minimum(map(p -> p.y,vert_func_vals_at_intersection))
    # end
    min_horz_func_val = minimum(map(p -> p.y,horz_func_vals_at_intersection))
    # min_vert_func_val = z_min+z_length
    vertfunccut = vert_cutout(vert_coord, z_min)
    
    ~,~,vert_func_vals_at_intersection = masked_occluded_intersections(vert_planks[j],vertfunccut)
    if (isempty(vert_func_vals_at_intersection))
      println("$j length of vert_plankj is $(length(vert_planks[j]))")
      

      #PointsNCurves.DrawSVGCurves("vert_curves_$j.svg",map(c -> Curve(c), [vert_planks[j]]))
    end

    min_vert_func_val = minimum(map(p -> p.y,vert_func_vals_at_intersection))
    
    cut_depth = depth_factor * min(min_horz_func_val, min_vert_func_val)
    
    if (cut_depth < z_min/2)
       error("somehow scaling the function values")
    end
      horz_cut = horz_cutout(horz_coord, cut_depth)
      push!(horz_cuts,horz_cut)
      (~, occluded_horz_plank, horz_intersections) = masked_occluded_intersections(horz_planks[i],horz_cut)
      if (length(occluded_horz_plank) <= 2)
        piece1 = occluded_horz_plank[1]
        if (length(occluded_horz_plank) == 2)
          piece2 = occluded_horz_plank[2]
          cut_piece = cleancurve(vcat(piece1,simple_cut(piece1[end],piece2[1],cut_depth),piece2))
          horz_planks[i] = ClosedCurve(cut_piece)
        else
          cut_piece = cleancurve(vcat(piece1,simple_cut(piece1[end],piece1[1],cut_depth)))
          horz_planks[i] = ClosedCurve(cut_piece)
        end
      else
        println(occluded_horz_plank[1])
        println("\n\n\n")
        println(occluded_horz_plank[2])
        error("occluded_horz_plank has length $(length(occluded_horz_plank))")
      end


      vert_cut = vert_cutout(vert_coord, cut_depth)
      push!(vert_cuts,vert_cut)
      (~, occluded_vert_plank, vert_intersections) = masked_occluded_intersections(vert_planks[j],vert_cut)
      if (length(occluded_vert_plank) <= 2)
        piece1 = occluded_vert_plank[1]
        if (length(occluded_vert_plank) == 2)
          piece2 = occluded_vert_plank[2]
          cut_piece = cleancurve(vcat(piece1,simple_cut(piece1[end],piece2[1],cut_depth),piece2))
          vert_planks[j] = ClosedCurve(cut_piece)
        else
          cut_piece = cleancurve(vcat(piece1,simple_cut(piece1[end],piece1[1],cut_depth)))
          vert_planks[j] = ClosedCurve(cut_piece)
        end
      else
        # glue_curve = []
        # for h in occluded_vert_plank
        #   # println("h has length $(length(h))")
        #   glue_curve = vcat(glue_curve, h)
        # end
        glue_curve = reduce(vcat, occluded_vert_plank)
        vert_planks[j] = ClosedCurve(glue_curve)
        DrawSVGCurves("vert_curves-1.svg",map(c -> Curve(c), [glue_curve]))
        DrawSVGCurves("vert_curves0.svg",map(c -> Curve(c), occluded_vert_plank))
        DrawSVGCurves("vert_curves1.svg",map(c -> Curve(c), vcat(occluded_vert_plank,vert_cut)))
        # for v in occluded_vert_plank
        #   println(v)
        #   println("\n\n\n") 
        # end
        # error("occluded_vert_plank has length $(length(occluded_vert_plank))")
      end
    end # inner for
  end # outer for
  (horz_planks,vert_planks)
end

function sample_grid_func(f::Function,Nvert::Int,Nhorz::Int,fsvert::Int, fshorz::Int,vert_min::Float64,vert_max::Float64,horz_min::Float64,horz_max::Float64)
 (Nvert < 3) && error("Need more vertical lines. the sample domain 
    should have 2 more than number of planks in a given direction")
  (Nhorz < 3) && error("Need more horizontal lines. the sample domain 
    should have 2 more than number of planks in a given direction")


  num_horz_sample_points = fshorz
  num_vert_sample_points = fsvert

  horz_x_coords = linspace(horz_min, horz_max, num_horz_sample_points)
  vert_x_coords = linspace(vert_min, vert_max, num_vert_sample_points)

  # "x" values that are constant for the vertical curves
  vert_cross_domain_vals = linspace(horz_min, horz_max, Nvert)[2:end-1]

  # "y" values that are constant for the horizontal curves
  horz_cross_domain_vals = linspace(vert_min, vert_max, Nhorz)[2:end-1]

  vert_curves = Curves()
  horz_curves = Curves()

  for (crv_idx, cross_val) in enumerate(vert_cross_domain_vals)
    crv = []
    for xval in vert_x_coords
      push!(crv,Point(xval, f(cross_val,xval)))
    end
    push!(vert_curves, crv)
  end

  for (crv_idx,cross_val) in enumerate(horz_cross_domain_vals)
    crv = []
    for xval in horz_x_coords
      push!(crv,Point(xval, f(xval,cross_val)))
    end
    push!(horz_curves, crv)
  end

  horz_curves, vert_curves
end

function glowify_grid_func(f::Function,sample_domain::DiscreteDomain,mat_width::Float64,horz_length::Float64,vert_length::Float64,z_length::Float64,z_min::Float64,fs::Int=100)
  (sample_domain.Nx < 3) && error("Need more vertical lines. the sample domain 
    should have 2 more than number of planks in a given direction")
  (sample_domain.Ny < 3) && error("Need more horizontal lines. the sample domain 
    should have 2 more than number of planks in a given direction")

  num_horz_sample_points = round.(Int,ceil(fs*horz_length))
  num_vert_sample_points = round.(Int,ceil(fs*vert_length))

  vert_min = sample_domain.minmin.y
  vert_max = sample_domain.maxmax.y
  horz_min = sample_domain.minmin.x
  horz_max = sample_domain.maxmax.x
  
  unscaled_horz_curves, unscaled_vert_curves = sample_grid_func(f,sample_domain.Nx, sample_domain.Ny, num_vert_sample_points, num_horz_sample_points, vert_min, vert_max,horz_min, horz_max)

  horz_x_coords = map(p -> p.x, unscaled_horz_curves[1]) # Assuming the domain is the same accross horizontal curves
  vert_x_coords = map(p -> p.x, unscaled_vert_curves[1]) # Assuming the domain is the same accross horizontal curves
  
  # The following values will be used to compute the sample points and cut points
  function rescale(vec,scl) 
    m = minimum(vec)
    M = maximum(vec)
    map(v -> scl * ((v - m)/(M - m)), vec)
  end
  scaled_horz_x_vals = rescale(horz_x_coords, horz_length)
  scaled_vert_x_vals = rescale(vert_x_coords, vert_length)

  vert_curves = Curves()
  horz_curves = Curves()

  for (i,c) in enumerate(unscaled_horz_curves)
    horz_func_vals = map(p -> p.y, c)
    crv = Vector{Point}()
    for (px,py) in zip(scaled_horz_x_vals, horz_func_vals)
      push!(crv, Point(px,py))
    end
    push!(horz_curves, crv)
  end

  for (i,c) in enumerate(unscaled_vert_curves)
    vert_func_vals = map(p -> p.y, c)
    crv = Vector{Point}()
    for (px,py) in zip(scaled_vert_x_vals, vert_func_vals)
      push!(crv, Point(px,py))
    end
    push!(vert_curves, crv)
  end
  
  glowify_curves_on_grid(vert_curves, horz_curves, mat_width, horz_length, vert_length, z_length, z_min)
end

#       5) translate the curves so that they don't overlap
#       6) export to svg

#function disjointify(planks)

#function export_to_svg(disjoint_planks,filename)

