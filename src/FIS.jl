# Contains FIS and Rule types
# ---------------------------
#
'''
	# A Rule connecting input and output membership function
	#
	# Properties
	# ----------
	# `input_mf_names` is a Vector of AbstractString (in order of inputs) which contains name of membership functions
	# `output_mf` is either an AbstractString containing the connected output membership function name with `input_mf_names`
	# 		or is a Vector of AbstractFloat containing consequence parameters in case of Sugeno fis
'''
struct Rule
	input_mf_names::Vector{AbstractString}
	output_mf::Union{Vector{Float64}, AbstractString}
	firing_method::AbstractString

	function Rule(input_mf_names,output_mf,firing_method="MIN")
		new(input_mf_names,output_mf,firing_method)
	end
end

'''
	# Fuzzy Inference System of Mamdani type
	#
	# Properties
	# ----------
	# `input_mfs_dicts` is a Vector of Dict (one for each input) which contains membership function
	# 		associated with its name
	# 		e.g. "small" => TriangularMF(1, 2, 3)
	# `output_mfs_dict` is a Dict containing membership function of output associated with names
	# 		e.g. "small" => GaussianMF(4, 2)
	# `rules` is a Vector of Rule that form the inference system
	# `output_domain` Discourse universe for the output
'''
struct FISMamdani
	input_mfs_dicts::Vector{Dict{AbstractString, MF}}
	output_mfs_dict::Dict{AbstractString,MF}
	rules::Vector{Rule}
	output_domain::Tuple{Real,Real}
end

'''
	# Fuzzy Inference System of Mamdani type
	#
	# Properties
	# ----------
	# `input_mfs_dicts` is a Vector of Dict (one for each input) which contains membership function
	# 		associated with its name
	# 		e.g. "small" => TriangularMF(1, 2, 3)
	# `rules` is a Vector of Rule (of sugeno type) that form the inference system
'''
struct FISSugeno
	input_mfs_dicts::Vector{Dict{AbstractString, MF}}
	rules::Vector{Rule}
end
