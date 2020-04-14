using CSV
using DataFrames
using LinearAlgebra

function build_data_dictionary()

    data_dictionary = Dict{AbstractString,Any}()

    stoichiometric_matrix = DataFrame(CSV.File("Network.csv",header=false))
    stoichiometric_matrix = convert(Matrix,stoichiometric_matrix)
    data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix

    E = 0.01E-3 # (Units: mmol/gDW)

    # Metabolic V-max for v1 -> v5
    v1_max = (203)*(3600)*(E) # v1 ATP + L-Citrulline + L-Aspartate --> AMP + Diphosphate + N-(L-Arginino)succinate (Units: 1/hr)
    v2_max = (34.5)*(3600)*E  # v2 N-(L-Arginino)succinate --> Fumarate + L-Arginine (Units: 1/hr)
    v3_max = (249)*(3600)*(E) # v3 L-Arginine + H2O --> L-Ornithine + Urea (Units: 1/hr)
    v4_max = (88.1)*(3600)*(E) # v4 Carbamoyl_phosphate + L-Ornithine --> Orthophosphate + L-Citrulline (Units: 1/hr)
    v5f_max = (13.7)*(3600)*(E) # v5f 2.0*L-Arginine + 4.0*Oxygen + 3.0*NADPH + 3.0*H --> 2.0*Nitric_oxide + 2.0*L-Citrulline + 3.0*NADP + 4.0*H2O (Units: 1/hr)
    v5r_max = (13.7)*(3600)*(E) # v5r 2.0*Nitric_oxide + 2.0*L-Citrulline + 3.0*NADP + 4.0*H2O --> 2.0*L-Arginine + 4.0*Oxygen + 3.0*NADPH + 3.0*H (Units: 1/hr)

    cell_mass = 2.3E-9 #cell mass of HeLa cell obtained from bionumber (Unit: g)
    cell_volume = 2.45*E-6 #cell volume of HeLa cell obtained from bionumber (Unit: m3)
    water = 0.7 #ratio of water in Hela cell obtained from bionumber in percentage (70%)
    cell_dry_mass = (1-water)*cell_mass
    cell_dry_weight = cell_volume/cell_dry_mass
    Conc_aspartate = 1.49E-2*cell_dry_weight #metabolite concentration of aspartate
    Conc_fumarate = 4.85E-4*cell_dry_weight #metabolite concentration of fumarate
    Conc_arginine = 2.55E-4*cell_dry_weight #metabolite concentration of arginine
    Conc_ornithine = 4.49E-3*cell_dry_weight #metabolite concentration of Ornithine (Saccharomyces cerevisiae)
    Conc_ATP = 4.67E-3*cell_dry_weight #metabolite concentration of ATP
    Conc_AMP = 4.23E-5*cell_dry_weight #metabolite concentration of AMP
    Conc_NADPH = 6.54E-5*cell_dry_weight #metabolite concentration of NADPH
    Conc_NADP = 2.84E-5*cell_dry_weight #metabolite concentration of NADP
    Conc_carbomyl_phosphate = 0 #metabolite concentration of carbomyl_phosphate is not available
    Conc_citruniline = 0 #metabolite concentration of citruniline is not available
    Conc_arginosuccinate = 0 #metabolite concentration of arginosuccinate is not available
    Conc_Urea = 0 #metabolite concentration of Urea is not available
    Conc_phosphate = 0 #metabolite concentration of phosphate is not available
    Conc_diphosphate = 0 #metabolite concentration of diphosphate is not available
    Conc_H2O = 0 #metabolite concentration of H2O is not available
    Conc_H = 0 #metabolite concentration of H is not available
    Conc_O2 = 0 #metabolite concentration of O2 is not available
    Conc_NO = 0 #metabolite concentration of NO is not available

    Km_aspartate = 1.54E-4*cell_dry_weight #metabolite concentration of aspartate
    Km_fumarate = 5.3E-3*cell_dry_weight #metabolite concentration of fumarate
    Km_arginine_v5 = 3E-3*cell_dry_weight #metabolite concentration of arginine for v2
    Km_arginine_v3 = 1.55E-3*cell_dry_weight #metabolite concentration of arginine for v3
    Km_ornithine = 1.6E-3*cell_dry_weight #metabolite concentration of Ornithine (Saccharomyces cerevisiae)
    Km_ATP = 3.92E-4*cell_dry_weight #metabolite concentration of ATP (Mus musculus)
    Km_AMP = 0 #metabolite concentration of AMP not available
    Km_NADPH = 0 #metabolite concentration of NADPH not available
    Km_NADP = 0 #metabolite concentration of NADP not available
    Km_carbomyl_phosphate = 0 #metabolite concentration of carbomyl_phosphate is not available
    Km_citruniline = 0 #metabolite concentration of citruniline is not available
    Km_arginosuccinate = 0 #metabolite concentration of arginosuccinate is not available
    Km_Urea = 0 #metabolite concentration of Urea is not available
    Km_phosphate = 0 #metabolite concentration of phosphate is not available
    Km_diphosphate = 0 #metabolite concentration of diphosphate is not available
    Km_H2O = 0 #metabolite concentration of H2O is not available
    Km_H = 0 #metabolite concentration of H is not available
    Km_O2 = 0 #metabolite concentration of O2 is not available
    Km_NO = 0 #metabolite concentration of NO is not available

    # Metabolic V for v1 -> v5
    v_1 = v1_max*((Conc_ATP)/(Conc_ATP + Km_ATP))*((Conc_aspartate)/(Conc_aspartate + Km_aspartate))
    v_2 = v2_max*1
    v_3 = v3_max*((Conc_arginine)/(Conc_arginine + Km_arginine_v3))
    v_4 = v4_max*((Conc_ornithine)/(Conc_ornithine+ Km_ornithine))
    v_5f = v5f_max*((Conc_arginine)/(Conc_arginine + Km_arginine_v5))
    v_5r = v5r_max*(1)

    # Building the Metabolic Flux Bounds Array
    metabolic_flux_bounds_array = [
    0.0 v_1;       # v1  (Units: mmol/gDW-hr)
    0.0 v_2;       # v2  (Units: mmol/gDW-hr)
    0.0 v_3;       # v3  (Units: mmol/gDW-hr)
    0.0 v_4;       # v4  (Units: mmol/gDW-hr)
    0.0 v_5f;      # v5f (Units: mmol/gDW-hr)
    -v_5r 0.0;     # v5r (Units: mmol/gDW-hr)
    0.0 10.0;      # b1 [] -> Carbamoyl_phosphate (Units: mmol/gDW-hr)
    0.0 10.0;      # b2 [] -> L-Aspartate (Units: mmol/gDW-hr)
    0.0 10.0;      # b3 Fumarate -> [] (Units: mmol/gDW-hr)
    0.0 10.0;      # b4 Urea -> [] (Units: mmol/gDW-hr)
    0.0 10.0;      # b5 [] -> ATP (Units: mmol/gDW-hr)
    0.0 10.0;      # b6 AMP -> [] (Units: mmol/gDW-hr)
    0.0 10.0;      # b7 Diphosphate -> [] (Units: mmol/gDW-hr)
    0.0 10.0;      # b8 Orthophosphate -> [] (Units: mmol/gDW-hr)
    0.0 10.0;      # b9 [] -> Oxygen (Units: mmol/gDW-hr)
    0.0 10.0;      # b10 [] -> NADPH (Units: mmol/gDW-hr)
    0.0 10.0;      # b11 [] -> H (Units: mmol/gDW-hr)
    0.0 10.0;      # b12 Nitric_oxide -> [] (Units: mmol/gDW-hr)
    0.0 10.0;      # b13 NADP -> [] (Units: mmol/gDW-hr)
    0.0 10.0;      # b14 [] -> H20 (Units: mmol/gDW-hr)
    0.0 10.0;      # b15 H20 -> [] (Units: mmol/gDW-hr)
    ]

    # Setup default flux bounds
    data_dictionary["metabolic_flux_bounds_array"] = metabolic_flux_bounds_array

    # Setup default species bounds array
    species_bounds_array = [
    0.0 0.0; #1 ATP
    0.0 0.0; #2 L-Citrulline
    0.0 0.0; #3 L-Aspartate
    0.0 0.0; #4 AMP
    0.0 0.0; #5 Diphosphate
    0.0 0.0; #6 N-(L-Arginino)succinate
    0.0 0.0; #7 Fumarate
    0.0 0.0; #8 L-Arginine
    0.0 0.0; #9 H20
    0.0 0.0; #10 L-Ornithine
    0.0 0.0; #11 Urea
    0.0 0.0; #12 Carbamoyl_phosphate
    0.0 0.0; #13 Orthophosphate
    0.0 0.0; #14 Oxygen
    0.0 0.0; #15 NADPH
    0.0 0.0; #16 H
    0.0 0.0; #17 Nitric_oxide
    0.0 0.0; #18 NADP
    ]

    data_dictionary["species_bounds_array"] = species_bounds_array

    # Setup the objective coefficient array
    objective_coefficient_array = [
    0.0;       # v1 (Units: mmol/gDW-hr)
    0.0;       # v2 (Units: mmol/gDW-hr)
    0.0;       # v3 (Units: mmol/gDW-hr)
    0.0;       # v4 (Units: mmol/gDW-hr)
    0.0;       # v5f (Units: mmol/gDW-hr)
    0.0;       # v5r (Units: mmol/gDW-hr)
    0.0;       # b1 [] -> Carbamoyl_phosphate (Units: mmol/gDW-hr)
    0.0;       # b2 [] -> L-Aspartate (Units: mmol/gDW-hr)
    0.0;       # b3 Fumarate -> [] (Units: mmol/gDW-hr)
    1.0;       # b4 Urea -> [] (Units: mmol/gDW-hr)
    0.0;       # b5 [] -> ATP (Units: mmol/gDW-hr)
    0.0;       # b6 AMP -> [] (Units: mmol/gDW-hr)
    0.0;       # b7 Diphosphate -> [] (Units: mmol/gDW-hr)
    0.0;       # b8 Orthophosphate -> [] (Units: mmol/gDW-hr)
    0.0;       # b9 [] -> Oxygen (Units: mmol/gDW-hr)
    0.0;       # b10 [] -> NADPH (Units: mmol/gDW-hr)
    0.0;       # b11 [] -> H (Units: mmol/gDW-hr)
    0.0;       # b12 Nitric_oxide -> [] (Units: mmol/gDW-hr)
    0.0;       # b13 NADP -> [] (Units: mmol/gDW-hr)
    0.0;       # b14 [] -> H20 (Units: mmol/gDW-hr)
    0.0;       # b15 H20 -> [] (Units: mmol/gDW-hr)
    ]

    data_dictionary["objective_coefficient_array"] = objective_coefficient_array

    data_dictionary["min_flag"] = false

    return data_dictionary

end
