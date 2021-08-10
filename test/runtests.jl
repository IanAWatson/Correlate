using Correlate
using Test
using Statistics

# We have run correlate and want to see if the result is within tolerance.

function check_correlation(zdata::Array{Float64,2}, corm::Array{Float64,2}, tol::Float64)
  (nrows,ncols) = size(zdata)

  cor_after = cor(zdata)

  for c in 1:ncols
    for r in (c+1):ncols
      if abs(corm[r,c] - cor_after[r,c]) > tol
        return false
      end
    end
  end

  return true
end

function test1()
  data = rand(1000,3)
  corm = [1 0.2 -0.3; 0.2 1 0.85; -0.3 .85 1]
  tol = 0.01
  rc = correlate!(data, corm, tol)
  if ! rc
    return false
  end

  if ! check_correlation(data, corm, tol)
    return false
  end

  (nrows,ncols) = size(data)

# now reverse it

  for c in 1:ncols
    for r in 1:ncols
      if r != c
        corm[r,c] = - corm[r,c]
      end
    end
  end

# write(STDERR, "CORM $(corm)\n")

  if ! correlate!(data, corm, tol)
    return false
  end

  if ! check_correlation(data, corm, tol)
    return false
  end

# now send back to uncorrelated

  for c in 1:ncols
    for r in 1:ncols
      if r != c
        corm[r,c] = 0.0
      end
    end
  end

  if ! correlate!(data, corm, tol)
    return false
  end

  return check_correlation(data, corm, tol)
end

# write your own tests here
@test test1()
