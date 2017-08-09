%BSMART Toolbox
%26-Feb-2007
%





%Aic test                  -  test Akaike Information Criterion(AIC) as a function of model order
%Whiteness test            -  test the noise process is a white(uncorrelated) process
%Consistency test          -  demonstrate the similarity between the correlation structure of the real data and that of AMVAR model simulated data
%Stability test            -  test the MVAR process is stable

%AMAR one window(pairwise) -  use AMAR method to analyse one window data(pairwise)
%AMAR moving window        -  use AMAR method to analyse data with moving window (pairwise)


%Granger Causality         -  caculate granger causality
%View data                 -  view original data by Chart View and Grid View
%Coherence View            -  view Coherence
%Granger Causality View    -  view granger causality
%Power View                -  view power
%Coherence network         -  view network of coherence
%Granger Causality network -  view network of granger causality
%Help
%   AIC_order_view             - M-file for AIC_order_view.fig
%   Lyapunov                   - M-file for Lyapunov.fig
%   MAR_coh3D                  - The only difference with coh3D is without global variables and t, f 
%   MAR_csd4dmem               - calculate spectral matrix, csd4d, based on output of MARfit 
%   MAR_make                   - To make it more general for any channel numbers
%   aic                        - M-file for aic.fig
%   aic_test                   - Compute the AIC
%   armorf                     - AR parameter estimation via LWR method by Morf modified.

%   autopower                  - Compute the auto power
%   bi_coherence               - Compute the pair coherence from the Bivariate models
%   bi_power                   - Compute the auto power from the Bivariate models
%   bsmart                     - MAR M-file for mar.fig

%   chartview                  - Chart view of data set
%   co_view                    - View coherence
%   coherence                  - M-file for coherence.fig
%   coherence_network          - M-file for coherence_network.fig
%   coherence_view             - M-file for coherence_view.fig
%   conetwork                  - Analysis and visualization of the coherence network
%   consistency_test           - M-file for consistency_test.fig
%   consistencytest            - Consistency test
%   constchk                   - function [R, Ru, Rd, trueR, denoRm] = consistncycheck1(fname, mu);
%   contents                   - BSMART Toolbox
%   dataview                   - M-file for dataview.fig
%   draw_layout                - Draws a layout for a graph

%   ensembding                 - Generate ensemble dat of size N 
%   fpeak                      - Author:    Geng Jun
%   ga_view                    - View Granger causality
%   ganetwork                  - Analysis and visualization of the Granger causality network



%   granger_causality_network  - M-file for granger_causality_network.fig
%   granger_causality_view     - M-file for granger_causality_view.fig
%   grid_view                  - M-file for grid_view.fig
%   gridview                   - Grid view of data set
%   gtm_dist                   - Calculate the squared distances between two sets of data points. 
%   iconRead                   - read an image file and convert it to CData for a HG icon.
%   lyap                       - calculating Lyapunov Exponent of MVAR
%   lyap_batch                 - A batch version to compute Lyapunov exponen
%   make_layout                - Creates a layout from an adjacency matrix

%   mar_fft                    - M-file for mar_fft.fig
%   mar_gen                    - [Y] = mar_gen (mar, T)
%   mar_init                   - To make it more general for any channel numbers
%   mar_spectra_new            - The only difference between current version and mar_spectra.m is
%   minor                      - calculate a minor of matrix, X, corresponding to the element X_ij
%   mov_bi_ga                  - Compute the granger causality from the moving window Bivariate models
%   mov_bi_model               - Moving window for bivariate models
%   mov_mul_model              - Moving window for multivariate model
%   moving_window              - M-file for moving_window.fig
%   moving_window_pairwise     - M-file for moving_window_pairwise.fig
%   movingwin                  - function [R1,R2,diffR,ratio]=movingwin(gelb, stim, shif, arcoeff, arnoise);
%   mul_coh                    - Squared Multiple Coherence(the ith chan and others) calculation 

%   one_bi_ga                  - Compute the granger causality from the one window Bivariate models
%   one_bi_model               - One window for bivariate model
%   one_mul_model              - One window for multivariate model
%   one_window                 - M-file for one_window.fig
%   one_window_pariwise        - M-file for one_window_pariwise.fig
%   pair_coh                   - ordinary coherence

%   paircoherence              - Compute the coherence pairs
%   pairwise_coherence         - M-file for pairwise_coherence.fig
%   pairwise_granger_causality - M-file for pairwise_granger_causality.fig
%   part_coh                   - Squared Partial Coherence(the ith and jth channel) calculation 
%   part_spect                 - Calculate partial auto/cross-spectr between channel i and j conditoned on 


%   po_view                    - View power
%   power_pairwise             - M-file for power_pairwise.fig
%   power_view                 - M-file for power_view.fig
%   pre_sube                   - Subtract the ensemble mean
%   pre_sube_divs              - Subtract the ensemble mean and divide by standard deviation
%   pre_subt                   - Subtract the temporal mean
%   pre_subt_divs              - Subtract the temporal mean and divide 
%   preprocessing              - M-file for preprocessing.fig
%   pwcausal                   - Using Geweke's method to compute the causality between any two channels
%   readdat                    - Read data in binary format and convert them into Matlab data
%   readdata                   - M-file for readdata.fig
%   spectrum                   - Get the coherence spectrum
%   spectrum_analysis          - M-file for spectrum_analysis.fig



%   textbox                    - Draws A Box around the text 
%   textoval                   - Draws an oval around text objects
%   toposort                   - A Topological ordering of nodes in a directed graph
%   whiteness                  - M-file for whiteness.fig
%   whiteness_test             - Whiteness test
%   writedat                   - Write data in format that steve required
%   writedata                  - M-file for writedata.fig
