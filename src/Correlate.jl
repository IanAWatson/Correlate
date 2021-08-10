__precompile__()

module Correlate

export correlate!

using Printf
using Statistics

# fmt is ignored, because the arguments to printf must be a literal string. 
# Leave the door open if that ever changes

function _print_with_space(output::Any, fmt::AbstractString, v::Float64)
  if v >= 0.0
    print(output, " ")
  end
  @printf output " %.4f" v
end

# Given a matrix of differences between the current correlation and the desired
# correlation, `cdiff`, and a matrix of tolerances, return the number of columns
# out of tolerance, as well as the largest difference.
function _count_not_converged(cdiff::Array{Float64,2}, tol::Array{Float64,2})
  (nrows,ncols) = size(cdiff)

  rc = 0;
  maxdiff = 0.0

  for c1 in 1:ncols
    for c2 in (c1+1):ncols
      d = cdiff[c1,c2]
      if d > tol[c1,c2]
        if d > maxdiff
          maxdiff = d
        end
        rc += 1
      end
    end
  end

  return (rc, maxdiff)
end

# Choose two distinct items in the range [1:n].
function choose_two(n::Int64)
  while true
    r1 = rand(1:n)
    r2 = rand(1:n)
    if r1 != r2
      return (r1, r2)
    end
  end
end

function choose_two_ordered(n::Int64)
  (r1,r2) = choose_two(n)
  if (r1 < r2)
    return (r1, r2)
  else
    return (r2, r1)
  end
end

# Assumes (but does not check) that the input DATA is zero mean within each column
# This is the function that does all the work

function correlate_zm!(data::Array{Float64, 2}, target_correlation::Array{Float64, 2},
                       tol::Array{Float64, 2}; maxiter = 0, report = 0, maxsec = 0,
                       verbose = false)

# write(stderr, "$(Base.corm(data, 0.0))\n")
  (nrows,ncols) = size(data)

  # Diagonal and cross terms of the correlation computation.
  xx = zeros(Float64, ncols)
  xy = zeros(Float64, ncols, ncols)

  for c in 1:ncols
#   write(stderr, "mean of column $(c) $(mean(data[:,c]))\n")
    for r in 1:nrows
      xx[c] += (data[r,c] * data[r,c])
    end
  end

  for c in 1:ncols
    xx[c] = sqrt(xx[c])
  end

  for c1 in 1:ncols
    for c2 in (c1+1):ncols
      for r in 1:nrows
        xy[c1,c2] += (data[r,c1] * data[r,c2])
      end
    end
  end

  for c1 in 1:ncols
    for c2 in (c1+1):ncols
      xy[c2,c1] = xy[c1,c2]
    end
  end

  mycor = zeros(Float64, ncols, ncols)

  for c1 in 1:ncols
    mycor[c1,c1] = 1.0
    for c2 in (c1+1):ncols
      mycor[c1,c2] = mycor[c2,c1] = xy[c1,c2] / (xx[c1] * xx[c2])
    end
  end

# write(stderr, "My correlation $(mycor)\n")
# write(stderr, "   correlation $(cor(data,mean=0))\n")

  cdiff = abs.(target_correlation - mycor)

  (nc,md) = _count_not_converged(cdiff, tol)

  if 0 == nc
    println(stderr, "Already converged")
    return true
  end

  iter = 0

  if 0 == maxiter
    maxiter = 100 * nrows
  end

  next_report = 0 == report ? (maxiter+maxiter) : report

  tzero = time()

  switches_done = 0

  switches_last_report = 0

  last_switch = 0

  column_switched = zeros(Int64, ncols)

  converged = false

  while iter < maxiter && ! converged
    iter += 1

    if iter > next_report
      (nc,md) = _count_not_converged(cdiff, tol)
      @printf(stderr, "%d iterations, %d switches (%d this period), %d not converged, max diff %f\n", iter, switches_done, (switches_done - switches_last_report), nc, md)
      if maxsec > 0 && time() - tzero > maxsec
        @printf(stderr, "Timeout after %f seconds\n", time() - tzero)
        break
      end
      next_report += report
      switches_last_report = switches_done
    end

    (row1,row2) = choose_two_ordered(nrows)

    col_to_switch = improves_correlations(data, row1, row2, target_correlation, xx, xy, cdiff, tol)

    if col_to_switch < 1 
      if iter - last_switch > div(maxiter, 10)
        @printf(stderr, "No progress in %d steps\n", iter - last_switch)
        break
      end
      continue
    end

    column_switched[col_to_switch] += 1

    last_switch = iter

    converged = true      # until we find otherwise

    for c1 in 1:ncols
      for c2 in 1:ncols
        if c1 == c2 continue end
        if c1 == col_to_switch || c2 == col_to_switch
          xy[c1,c2] += - data[row1,c1]*data[row1,c2] - data[row2,c1]*data[row2,c2] + 
                       data[row1,c1]*data[row2,c2] + data[row1,c2]*data[row2,c1]
          mycor[c1,c2] = xy[c1,c2]/(xx[c1] * xx[c2])
          cdiff[c1,c2] = abs(target_correlation[c1,c2] - mycor[c1,c2])
          if converged && cdiff[c1,c2] > tol[c1,c2]
            converged = false
          end
        elseif cdiff[c1,c2] > tol[c1,c2]
          converged = false
        end
      end
    end

    data[row1,col_to_switch],data[row2,col_to_switch] = data[row2,col_to_switch],data[row1,col_to_switch]
    switches_done += 1
  end

  if verbose
    if converged 
      println(stderr, "Converged. ")
#     write(stderr, "cdiff $(cdiff)")
    end
    println(stderr, "Performed $iter iterations, made $switches_done switches")
    sum_errors = 0.0
    not_converged = 0
    for c1 in 1:ncols
      for c2 in 1:ncols
        _print_with_space(stderr, " %3f", mycor[c1,c2])
        if c2 > c1
          sum_errors += abs(target_correlation[c1,c2] - mycor[c1,c2])
          if abs(target_correlation[c1,c2] - mycor[c1,c2]) > tol[c1,c2]
#           @printf(stderr, "col %d vs col %d target %f mycor %f not converged\n", c1, c2, target_correlation[c1,c2], mycor[c1,c2])
            not_converged += 1
          end
        end
      end
      @printf(stderr, "     ")
      for c2 in 1:ncols
        _print_with_space(stderr, " %3f", target_correlation[c1,c2])
      end
      println(stderr)
    end
    println(stderr, "Sum of errors $sum_errors, $not_converged items not converged")
    for c in 1:ncols
      println(stderr, "Switched column $c $(column_switched[c]) times")
    end
    @printf(stderr, "Calculation took %.2f seconds\n", time()-tzero)
  end

  return converged
end

function improves_correlations(data::Array{Float64,2}, row1::Int64, row2::Int64, target_correlation::Array{Float64, 2},
                               xx::Array{Float64, 1}, xy::Array{Float64,2}, cdiff::Array{Float64,2},
                               tol::Array{Float64,2})
  (nrows,ncols) = size(data)
 
  current_diff = 0.0
 
  for c1 in 1:ncols
    for c2 in (c1+1):ncols
      d = cdiff[c1,c2]
      if d > tol[c1,c2]
        current_diff += d
      end
    end
  end
# @printf(stderr, "Rows %d and %d, cirrent diff %f\n", row1, row2, current_diff)

  lowest = current_diff
  col_to_switch = -1
 
  for c in 1:ncols     # what happens if we swap column C
    sum = 0.0
    for x1 = 1:ncols      # check all correlations involving C
      for x2 = (x1+1):ncols
        if x1 == c || x2 == c
          t = xy[x1,x2] - data[row1,x1]*data[row1,x2] - data[row2,x1]*data[row2,x2] + 
                          data[row1,x2]*data[row2,x1] + data[row1,x1]*data[row2,x2]
          tc = t /(xx[x1] * xx[x2])
#         jj =  xy[x1,x2]/(xx[x1] * xx[x2])

#         @printf(stderr, "switching col %d, testing %d and %d, c %f current %f target %f\n", c, x1, x2, tc, jj, target_correlation[x1,x2])
          d = abs(tc - target_correlation[x1,x2])
          if d > tol[x1,x2]
            sum += d
          end
        else
          d = cdiff[x1,x2]
          if d > tol[x1,x2]
            sum += d
          end
        end
      end
    end
#   @printf(stderr, "switching column %d, sum %f\n", c, sum)
    if sum < lowest
      lowest = sum
      col_to_switch = c
    end
  end
 
# @printf(stderr, "Existing diff %f, new %f, column %d\n", current_diff, lowest, col_to_switch)

  return col_to_switch
end

function improves_correlations(data::Array{Float64,2}, row1::Int64, row2::Int64, target_correlation::Array{Float64, 2},
                               xx::Array{Float64, 1}, xy::Array{Float64,2}, cdiff::Array{Float64,2})
  (nrows,ncols) = size(data)
 
  current_diff = 0.0
 
  for c1 in 1:ncols
    for c2 in (c1+1):ncols
      current_diff += cdiff[c1,c2]
    end
  end
# @printf(stderr, "Rows %d and %d, cirrent diff %f\n", row1, row2, current_diff)

  lowest = current_diff
  col_to_switch = -1
 
  for c in 1:ncols     # what happens if we swap column C
    sum = 0.0
    for x1 = 1:(ncols-1)      # check all correlations involving C
      for x2 = (x1+1):ncols
        if x1 == c || x2 == c
          t = xy[x1,x2] - data[row1,x1]*data[row1,x2] - data[row2,x1]*data[row2,x2] + 
                          data[row1,x2]*data[row2,x1] + data[row1,x1]*data[row2,x2]
          tc = t /(xx[x1] * xx[x2])
#         jj =  xy[x1,x2]/(xx[x1] * xx[x2])

#         @printf(stderr, "switching col %d, testing %d and %d, c %f current %f target %f\n", c, x1, x2, tc, jj, target_correlation[x1,x2])
          sum += abs(tc - target_correlation[x1,x2])
        else
          sum += cdiff[x1,x2]
        end
      end
    end
#   @printf(stderr, "switching column %d, sum %f\n", c, sum)
    if sum < lowest
      lowest = sum
      col_to_switch = c
    end
  end
 
# @printf(stderr, "Existing diff %f, new %f, column %d\n", current_diff, lowest, col_to_switch)

  return col_to_switch
end

# Basically convert the data to zero mean, call the underlying function, then shift back

"""
    correlate!(data::Array{Float64,2}, corm::Array{Float64,2}, tolm::Array{Float64,2};
               maxiter=100*nrows, report=0, maxsec=0, verbose=false)

Enforce a correlation structure on columns of independently derived data.\n
CORM: matrix of desired correlations between columns.
TOLM: tolerance associated with each pair.

# Arguments

* maxiter: number of iterations (def 100*nrows).
* report: report progress every <report> iterations.
* maxsec: give up after <maxsec> seconds of wall clock time.
* verbose: write computation summary to stderr.

# Examples

    data = rand(1000,4);
    corm = fill(0.8, (4,4))     # diagonals should be 1, but are not checked
    tolm = fill(0.01, (4,4))    # diagonals never examined
    correlate!(data, corm, tolm, verbose=true)
    cor(data)      # check result, should show somewhere around 0.8

A variant where tolerance is the same for all correlations

    correlate!(data::Array{Float64,2}, corm::Array{Float64,2}, tol_value::Float64;
               maxiter=100*nrows, report=0, maxsec=0, verbose=false)

A variant for two columns, where correlation and tolerance are scalar

    correlate!(data::Array{Float64,2}, cor_value::Float64, tol_value::Float64;
               maxiter=100*nrows, report=0, maxsec=0, verbose=false)

Run chol(corm) on your proposed correlation matrix (must be unit diagnonals) to see if
a solution is feasible. Correlate! does not check this, but will yield a "partial"
result if given an impossible task.
"""
function correlate!(data::Array{Float64, 2}, corm::Array{Float64, 2}, tolm::Array{Float64, 2};
                    maxiter = 0, report = 0, maxsec = 0, verbose = false)
  (nrows,ncols) = size(data)

  column_mean = [mean(data[:,c]) for c in 1:ncols]

  for c in 1:ncols
    cmean = column_mean[c]
    for r in 1:nrows
      data[r,c] -= cmean
    end
  end

  rc = correlate_zm!(data, corm, tolm, maxiter=maxiter, report=report, maxsec=maxsec, verbose = verbose)

  for c in 1:ncols
    cmean = column_mean[c]
    for r in 1:nrows
      data[r,c] += cmean
    end
  end

  return rc
end

# Version with two columns, so just one value for correlation and tolerance

function correlate!(data::Array{Float64, 2}, cor_value::Float64, tol_value::Float64;
                    maxiter = 0, report = 0, maxsec = 0, verbose = false)
  (nrows,ncols) = size(data)
  if 2 != ncols
    @printf(stderr, "Correlate.correlate:version with scalar correlation and tolerance operates only on N*2 data\n")
    return false
  end

  corm = fill(cor_value, (2,2))
  tolm = fill(tol_value, (2,2))

  return correlate!(data, corm, tolm, maxiter=maxiter, report=report,maxsec=maxsec, verbose = verbose)
end

function correlate!(data::Array{Float64, 2}, corm::Array{Float64,2}, tol_value::Float64;
                   maxiter = 0, report = 0, maxsec = 0, verbose = false)
  (nrows,ncols) = size(data)

  tolm = fill(tol_value, (nrows,ncols))

  return correlate!(data, corm, tolm, maxiter=maxiter, report=report, maxsec=maxsec, verbose = verbose)
end

end # module Correlate
