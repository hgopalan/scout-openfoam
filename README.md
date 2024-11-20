Urban microclimate solver in OpenFOAM with some engineering hacks: 

1. Handles steady and unsteady solvers with temperature and humidity within the same solver.
2. Porous media model for trees. Plant a million trees with ease by converting GIS point shape files to Openfoam compatible dictonary files. 
3. Modified fvdom with leaf area density model for trees.
4. Leaf energy balance
5. 1-D conduction model and surface energy balance for temperature.
6. Air-con, heat sources and water sprinkling can be included as source terms on wall boundary conditions. 
7. One-equation Wray-Agarwal model for all stratifications and ASL/ABL support. Homogenity is ensured for fetch > 50 km within 0,1 m/s.
8. Thermal comfort and WBGT outputs.
9. Converges for very low wind speeds < 0.2 m/s and clear-sky radiation.
10. Run UHI simulations by converting netCDF output to the Openfoam compatible dictonary files. 
