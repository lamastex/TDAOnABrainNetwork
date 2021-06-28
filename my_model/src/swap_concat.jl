using CSV, DataFrames, ProgressMeter, Glob, DelimitedFiles, DelimitedFiles, NPZ
import CSV

m_type = [
    "L1_DAC", "L1_DLAC", "L1_HAC", "L1_NGC-DA", "L1_NGC-SA", "L1_SLAC", "L23_BP",
    "L23_BTC", "L23_ChC", "L23_DBC", "L23_LBC", "L23_MC", "L23_NBC", "L23_NGC",
    "L23_PC", "L23_SBC", "L4_BP", "L4_BTC", "L4_ChC", "L4_DBC", "L4_LBC",
    "L4_MC", "L4_NBC", "L4_NGC", "L4_PC", "L4_SBC", "L4_SP", "L4_SS", "L5_BP",
    "L5_BTC", "L5_ChC", "L5_DBC", "L5_LBC", "L5_MC", "L5_NBC", "L5_NGC", "L5_SBC",
    "L5_STPC", "L5_TTPC1", "L5_TTPC2", "L5_UTPC", "L6_BP", "L6_BPC", "L6_BTC",
    "L6_ChC", "L6_DBC", "L6_IPC", "L6_LBC", "L6_MC", "L6_NBC", "L6_NGC", "L6_SBC",
    "L6_TPC_L1", "L6_TPC_L4", "L6_UTPC"
    ]


function morph_sizes(a)
  sizes = similar(a, Int64)

  directory = "../output/GB/general_reconstruction/"
  for (i, name) in enumerate(a)
    file = "$(directory)$(name)$(name).csv"
    sizes[i] = length(readlines(file))
  end
  return sizes
end

function concat_files(a)
  directory = "../output/GB/general_reconstruction/"

  sizes = morph_sizes(a)
  total = sum(sizes)
  A = zeros(Bool, total, total)

  @showprogress for i in (1:length(a))
    for j in (1:length(a))
      k_start = sum(sizes[1:i-1])
      l_start = sum(sizes[1:j-1])

      file = "$(directory)$(a[i])$(a[j]).csv"
      open(file) do f
        for (k, row) in enumerate(eachline(f))
          for (l, v) in enumerate(row[1:2:end])
            A[k_start + k, l_start + l] = v == '1'
          end
        end
      end

    end
  end
  npzwrite("../output/GB/gen_biol.npy", A)
end

concat_files(m_type)
