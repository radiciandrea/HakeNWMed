# HakeNWMed
Repository to store code from Radici et al., "Assessing fish-fishery dynamics from a spatially explicit metapopulation perspective reveals winners and losers in fisheries management" (unp.)

All codes stored here need MATLAB to be run

modello_demografico_8_GSA9.m is the calibrated metapopulation model which calculates stock abundances for each age class in each cell at each time step (1 year) and outputs (SSB = Spauning stock biomass; L = landings)

it loads the following data from the data repository:

n_celle_regioni_3 (table associating cells id and regions name);
celle_M_4_3_o (ordered vector of cells id);
Mc_4_3_o_perc (row normalized connectivity matrix;
R10_EFF9_var (coefficients describing future recruitment on GSA 10 and fishing effort trend in GSA 9);
persistenza_n (coefficients describing persistency of the nursery in every cell);
persistenza_sg (coefficients describing persistency of the spowing grounds in every cell);
area_celle_o (size of cells in km2);

And the following scripts:

coefficienti_fisheries.m (contains fishing coefficients);
plot_modello_demografico_8_GSA_9.m (plot the outputs of the model);

These last scripts load other data:

D9 (matrix containing the euclidean distance among cells);
cella_temp (the id of the cell);
fisheries_GSA_9 (calibrated coefficiencies related to fishery in the GSA 9 by metier)
