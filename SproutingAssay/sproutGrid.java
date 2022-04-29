/*
GRIFFIN LAB
LAUREN PRUETT, ALEX TAING
UNIVERSITY OF VIRGINIA
2D Sprouting Angiogenesis ABM - Publication 2022
*/

package SproutingAssay;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Rand;
import HAL.Util;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.Timestamp;
import java.util.ArrayList;

public class sproutGrid extends AgentGrid2D<sproutAgent> {


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////      PARAMETERS      /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // BATCH RUNS
    public final static boolean BATCH_RUN = true; // sets whether a single trial run is performed or a batch run.  Only BATCH RUNS use trials.
                                                  // NOTE: if single trial run, will only run on the first percentage in HEPARIN_PERCENTAGES
    public final static boolean EXPORT_DATA = true;
    public final static boolean EXPORT_TIME_DATA = true; // (note that EXPORT_DATA must also be true as well to export time data)
    public final static boolean EXPORT_HEAD_CELL_DISTANCE_DATA = true; // (note that EXPORT_DATA must also be true as well to export distance data)
    public final static int TRIALS = 1;  // number of trials (only for BATCH RUNS)
    public final static double[] HEPARIN_PERCENTAGES = {0.1};  // percent HEPARIN ISLANDS that will be inilialized in the scaffold
    public final static double FOLD_CHANGE_SAMPLE_TIME = 0.25 ; // sampling interval for fold change data export (HOURS)


    // VESSEL PARAMETERS FIXED!
    public static final int CULTURE_RADIUS_MICRONS = 140; // radius of the initial assay spheroid (MICRONS)
    public final static int SIGHT_RADIUS_MICRONS = 20; // how far a EC can sense (MICRONS)
    public final static double PERSISTENCY_TIME_HOURS = 3; // time between changes direction for EC migration (HOURS)
    public static final int MAX_ELONGATION_LENGTH_MICRONS = 40; // max elongation length of a EC (MICRONS)
    public final static int MIGRATION_RATE_MICRONS_PER_HOUR = 30; // rate of elongation for EC head cells (MICRONS/HOUR)
    public final static double BRANCH_DELAY_TIME_HOURS = 6; // required delay time between EC branching in (HOURS)

    // VESSEL PARAMETERS FOR SENSITIVITY ANALYSIS
    public final static double VEGF_SENSITIVITY = 0.025; // minimal VEGF required to be sensed by EC (0.03 baseline minimum VEGF to attract cell growth) (RELATIVE UNITS)
    public static double VESSEL_VEGF_INTAKE = 0.001; // amount of the current VEGF that is consumed when a blood vessel is present in a grid cell (RELATIVE UNITS)
    public final static double VEGF_DEGRADATION_RATE = 0.5; // Fraction of VEGF that is degraded every 90 minutes in accordance with VEGFs half life
    public final static double REQUIRED_VEGF_GRADIENT_DIFFERENCE = 0.005; // Difference of VEGF concentration required between a HEAD cell's location and its next division location in order to elongate (RELLATIVE UNITS)

    // BRANCHING PROBABILITY AND THRESHOLDS_ PROBABILITIES
    public final static double INITIAL_PERCENT_HEAD_CELLS = 0.07; // probability during spheroid initiation that any edge cell of the spheroid will be of type HEAD
    public final static double LOW_BRANCHING_PROBABILITY= 0.4; // probability of branching while VEGF is under LOW_MED_VEGF_THRESHOLD
    public final static double LOW_MED_VEGF_THRESHOLD = 0.05;  // the threshold between LOW and MED branching probability (RELATIVE UNITS)
    public final static double MED_BRANCHING_PROBABILITY= 0.55; // probability of branching while VEGF is between LOW_MED_VEGF_THRESHOLD and MED_HIGH_VEGF_THRESHOLD
    public final static double MED_HIGH_VEGF_THRESHOLD = 0.25; // the threshold between MED and HIGH branching probability (RELATIVE UNITS)
    public final static double HIGH_BRANCHING_PROBABILITY= 0.9; // probability of branching while VEGF is above MED_HIGH_VEGF_THRESHOLD

    // MAP GEL PARAMETERS - FIXED
    public final static int MAP_RADIUS_MICRONS = 40; // radius of MAP particles (MICRONS)
    public final static double MAP_SPACING_MICRONS = 15; // minimum space between the edges of two MAP particles (MICRONS)
    public final static double HEP_MAP_VEGF_RELEASE = 1; // amount of VEGF added per media exchange (RELATIVE UNITS) (NOTE: This is based on percentage of hep particles - HEP% * HEP_MAP_VEGF_RELEASE = 0.1)
    public final static double MEDIA_EXCHANGE_SCHEDULE_HOURS = 1; // exchange media to refresh VEGF every X hours

    // MAIN METHOD PARAMETERS - FIXED
    public final static int x_MICRONS = 4000; // x dimension of the wound-space (MICRONS)
    public final static int y_MICRONS = 4000; // y dimension of the wound-space (MICRONS)
    public final static int SCALE_FACTOR = 2; // for visualization: changes the scale of the pixels
    public final static int TICK_PAUSE = 1; // for model runtime: changes the amount of time between each tick
    public final static int RUNTIME_HOURS = 24; // the timeframe that the simulation models (HOURS)

    // DIFFUSION
    public final static double DIFFUSION_COEFFICIENT = 0.733; // diffusion coefficient (UNITLESS)



/////////////////////////////////////////////      CONVERSIONS      ////////////////////////////////////////////////////

    // SCALING FACTORS
    public final static int MICRONS_PER_PIXEL = 10; // 1 pixel represents 10 microns
    public final static int TICKS_PER_HOUR = 60; // 1 tick represents 1 minute

    // vessel unit conversions
    public static final int CULTURE_RADIUS = CULTURE_RADIUS_MICRONS/ MICRONS_PER_PIXEL;
    public final static int SIGHT_RADIUS = SIGHT_RADIUS_MICRONS/ MICRONS_PER_PIXEL;
    public static final int MAX_ELONGATION_LENGTH = MAX_ELONGATION_LENGTH_MICRONS/ MICRONS_PER_PIXEL;
    public final static double MIGRATION_RATE = 1/((MIGRATION_RATE_MICRONS_PER_HOUR/(double) MICRONS_PER_PIXEL)*(1/(double)TICKS_PER_HOUR)); // convert to "elongate every ___ ticks"
    public final static double PERSISTENCY_TIME = PERSISTENCY_TIME_HOURS * TICKS_PER_HOUR;
    public final static int BRANCH_DELAY_TIME = (int)(BRANCH_DELAY_TIME_HOURS * TICKS_PER_HOUR);

    // particle unit conversions
    public final static int MAP_RADIUS = MAP_RADIUS_MICRONS/ MICRONS_PER_PIXEL;
    public final static double MAP_SPACING = (double)(MAP_RADIUS_MICRONS)/ MICRONS_PER_PIXEL + MAP_SPACING_MICRONS/ MICRONS_PER_PIXEL;
    public final static double MEDIA_EXCHANGE_SCHEDULE_TICKS = MEDIA_EXCHANGE_SCHEDULE_HOURS*TICKS_PER_HOUR;

    // grid
    public final static int x = x_MICRONS/ MICRONS_PER_PIXEL;
    public final static int y = y_MICRONS/ MICRONS_PER_PIXEL;

    // runtime
    public final static int TIMESTEPS = RUNTIME_HOURS*TICKS_PER_HOUR;
    public final static int FOLD_CHANGE_SAMPLE_TICKS = (int)(FOLD_CHANGE_SAMPLE_TIME*TICKS_PER_HOUR); // take a sample every ____ ticks



//////////////////////////////////////////////      MISC VARIABLES      ////////////////////////////////////////////////

    public static int HEAD_CELL = sproutAgent.HEAD_CELL;  // Head cell type
    public static int BODY_CELL = sproutAgent.BODY_CELL;  // Body cell type

    Rand rng = new Rand();
    int[] divHood = Util.VonNeumannHood(false); // neighborhood for division: cells can only divide/elongate to adjacent cells (cardinal directions)
    int[] VEGFHood = Util.CircleHood(false, SIGHT_RADIUS); // the range that cells can see VEGF concentrations (for chemotaxis)
    int[] MAP_rad = Util.CircleHood(true, MAP_RADIUS); // radius of MAP particles
    int[] MAP_space = Util.CircleHood(true, MAP_SPACING); // "cushion" between MAP particles (dictates the the area around the center of a MAP particle of
                                                                     // of radius equal to MAP_RADIUS + MAP_SPACING)

    PDEGrid2D VEGF; // Initialize PDE Grid



////////////////////////////////////////////////      DATA EXPORT      /////////////////////////////////////////////////

    public static StringBuilder CSV = new StringBuilder();
    public static StringBuilder TIME_CSV = new StringBuilder();
    public static StringBuilder HEAD_CSV = new StringBuilder();
    public int initialCultureSize;
    public ArrayList<Double> foldChangeOverTime = new ArrayList<>();
    public ArrayList<Double> FCTime = new ArrayList<>();



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////      METHODS      ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    /**
     * Heparin Grid constructor. The grid created will be the vessel where all agents will reside.
     * @param x x dimension of grid
     * @param y y dimension of grid
     */
    public sproutGrid(int x, int y) { // Constructor for the agent grid
        super(x, y, sproutAgent.class);
        VEGF = new PDEGrid2D(x, y);
    }


    /**
     * Calls each cell to take actions by calling the StepCell() method of the sproutAgent class.
     * The sproutAgent class will then dictate what actions each cell will take depending on cell
     * type and its environment.
     *
     * The tick timer is incremented to advance the time.
     *
     * This function also calls the VEGF diffusion to iterate and update.  (without this, no diffusion
     * would occur.)
     */
    public void StepCells(){ // steps all the cells
        for (sproutAgent endoCell : this) {
            endoCell.StepCell();
        }
        IncTick();
        VEGF.DiffusionADI(DIFFUSION_COEFFICIENT); //Alternating Direction Implicit Method Diffusion
        VEGF.Update();
    }


    /**
     * Draws the PDE window.  This function provides purely visual functionality
     * by visualizing PDE diffusion of VEGF
     * @param windows the window that the visualization will be drawn in
     */
    public void DrawPDE(GridWindow windows){
        for (int i = 0; i < length; i++) {
            windows.SetPix(i,Util.HeatMapBGR(VEGF.Get(i)));
        }
    }


    /**
     * Draws the cell model.  This function provides visual functionality
     * by visualizing the agents (all cells and particles).
     * @param win the window that the visualization will be drawn in
     */
    public void DrawModel(GridWindow win){
        for (int i = 0; i < length; i++) {
            int color = Util.BLACK;
            if (GetAgent(i) != null) {
                sproutAgent cell = GetAgent(i);
                color = cell.color;
            }
            win.SetPix(i, color);
        }
    }


    /**
     * Replaces native Dist function which DOES NOT WORK with initVesselsCircleCulture
     * @param x1 first x coord
     * @param y1 first y coord
     * @param x2 second x coord
     * @param y2 second y coord
     * @return the integer distance between them
     */
    public int distance(double x1, double y1, double x2, double y2){
        double dist = Math.sqrt((Math.pow((x1-x2), 2)+Math.pow((y1-y2), 2)));
        return (int)dist;
    }


    /**
     * Initializes the spheroid for the assay
     * @param model the model to initialize the vessels in
     */
    public void initVesselsCircleCulture(sproutGrid model) {

        int center = I((Xdim()/2), (Ydim()/2));  // determines the center, where the spheroid will be placed
        for (int i = 0; i < (model.Xdim()*model.Ydim()); i++) {  // iterates through all coordinates
            double dist = distance(ItoX(i), ItoY(i), ItoX(center), ItoY(center));
            if (dist < (CULTURE_RADIUS)) { // if the distance for a coordinate is within CULTURE_RADIUS of the center
                if (dist > (CULTURE_RADIUS - 2)){  // and if the cell is on the edge of the spheroid, then
                    if (Math.random() < INITIAL_PERCENT_HEAD_CELLS) {
                        model.NewAgentSQ(i).InitVessel(sproutAgent.HEAD_CELL, 0);
                    }
                } else if (dist > (CULTURE_RADIUS-3)){
                    if (Math.random() < INITIAL_PERCENT_HEAD_CELLS) { // the cell may be initialized as type HEAD or BODY
                        model.NewAgentSQ(i).InitVessel(sproutAgent.HEAD_CELL, 0);
                    } else {
                        model.NewAgentSQ(i).InitVessel(sproutAgent.BODY_CELL, 0);
                    }
                } else { // else (if it is in within CULTURE_RADIUS but not on the edge), then
                    model.NewAgentSQ(i).InitVessel(sproutAgent.BODY_CELL, 0); // the cell will be of type BODY
                }
            }
        }
    }


    /**
     * Initializes MAP particles and Heparin microislands into the model in HCP pattern with
     * proper spacing as dictated by MAP_space
     * @param model model to initialize the particles in
     */
    public void initMAPParticles(sproutGrid model, double Heparin_Percent){

        for (int i = 0; i < x*y; i++) { // Iterate through every coordinate in the grid
            int cellType = sproutAgent.MAP_PARTICLE; // assume that it will be a MAP particle
            double chance = Math.random();
            if (chance < Heparin_Percent){ // if chosen probability dictates that it will he a heparin microIsland
                cellType = sproutAgent.HEPARIN_MAP;// then its type will be changed to heparin microIsland
            }

            int occlusions = MapOccupiedHood(MAP_space, i); // for the coordinate selected, check for occlusions (other cells) in area of radius MAP_space around it
            int open = MapEmptyHood(MAP_rad, i);
            if (occlusions == 0) { // if there are no occlusions
                for (int j = 0; j < open; j++){ // then initialize the particle at that location (checking that all particles stay within the window)
                    if (0 < MAP_rad[j] && MAP_rad[j] < x*y){
                        model.NewAgentSQ(MAP_rad[j]).Init(cellType);
                    }
                }
            }
        }
    }


    /**
     * Counts the total number of pixels which have a vessel present
     * (used in calculating fold change and total vessel length)
     * @return the total number of pixels that are occupied by vessels
     */
    public int countVessels() {
        int vessel_unit_counter = 0;
        for (int x_coord = 0; x_coord < x; x_coord++) {  // iterates through all coordinates
            for (int y_coord = 0; y_coord < y; y_coord++) {
                Iterable<sproutAgent> agents = IterAgents(x_coord, y_coord);
                for (sproutAgent agent : agents) {
                    if (agent.type == HEAD_CELL || agent.type == BODY_CELL){ // if one of the coordinates is a vessel
                        vessel_unit_counter ++; // tally and continue
                        break;
                    }
                }
            }
        }
        return vessel_unit_counter;
    }
    
    
    /**
     * Collects data stored in a static variable, CSV.
     * ExportData uses this collected information to export to a CSV file.
     */
    public void CollectData(double heparinPercentage){
        // vessel cells, MAP Particle, and Heparin MAP Data
        CSV.append("\n");
        
        // Note percentage
        CSV.append((int)Math.round(heparinPercentage*100)).append("%,");
        
        // Total vessel length
        int numVessels = countVessels();
        CSV.append((numVessels-initialCultureSize)* MICRONS_PER_PIXEL).append(","); // total length of vessels
        
        // Fold change
        double foldChange = numVessels/(double)initialCultureSize;
        CSV.append(foldChange).append(",");

        // Head cell data
        ArrayList<Double> headCellDistances = FinalHeadCellDistance();

        //grid size
        CSV.append(x_MICRONS).append(",");

        // diffusion coefficient
        CSV.append(DIFFUSION_COEFFICIENT).append(",");

        //vegf sensitivity
        CSV.append(VEGF_SENSITIVITY).append(",");

        //initial percent head cells
        CSV.append(INITIAL_PERCENT_HEAD_CELLS).append(",");

        //vessel vegf intake
        CSV.append(VESSEL_VEGF_INTAKE).append(",");

        //hours
        CSV.append(RUNTIME_HOURS).append(",");

        //vegf degradation
        CSV.append(VEGF_DEGRADATION_RATE).append(",");

        //Low-med threshold
        CSV.append(LOW_MED_VEGF_THRESHOLD).append(",");

        //med-high threshold
        CSV.append(MED_HIGH_VEGF_THRESHOLD).append(",");

        //VEGF in each hep particle
        CSV.append(HEP_MAP_VEGF_RELEASE).append(",");

        //How frequently hep-map vegf VEGF is released
        CSV.append(MEDIA_EXCHANGE_SCHEDULE_HOURS);


        // TIME DATA
        if (EXPORT_TIME_DATA){
            if (TIME_CSV.length() == 0){
                TIME_CSV.append("Time (hours), ");
                for (Double time : FCTime) {
                    TIME_CSV.append(time).append(",");
                }
                TIME_CSV.append("\n");
            }
            TIME_CSV.append("[").append((int) (heparinPercentage * 100)).append("%] Fold Change, ");
            for (Double TimeFoldChange : foldChangeOverTime) {
                TIME_CSV.append(TimeFoldChange).append(",");
            }
            TIME_CSV.append("\n");
        }

        // HEAD CELL DISTANCE DATA
        if (EXPORT_HEAD_CELL_DISTANCE_DATA){
            if (HEAD_CSV.length() == 0){
                HEAD_CSV.append("Final Head Cell Distances (microns) ").append("\n");
            }
            HEAD_CSV.append("[").append((int) (heparinPercentage * 100)).append("%],");
            for (Double distance : headCellDistances) {
                HEAD_CSV.append(distance).append(",");
            }
            HEAD_CSV.append("\n");
        }
    }


    /**
     * Calculates final head cell distances from the center of the wound.
     * @return an array of distances of each head cell to the center of the model
     */
    public ArrayList<Double> FinalHeadCellDistance(){
        // Head cell distances
        ArrayList<Double> head_cell_distances = new ArrayList<>();
        int center_x = (Xdim()/2);
        int center_y = (Ydim()/2);
        for (int x_coord = 0; x_coord < x; x_coord++) {
            for (int y_coord = 0; y_coord < y; y_coord++) {
                Iterable<sproutAgent> agents = IterAgents(x_coord, y_coord);
                for (sproutAgent agent : agents) {
                    if (agent.type == HEAD_CELL){
                        double dist = (distance(x_coord, y_coord, center_x, center_y))*MICRONS_PER_PIXEL;
                        head_cell_distances.add(dist);
                    }
                }
            }
        }

        return head_cell_distances;
    }


    /**
     * Exports data to a file labeled with run and date information to directory "SproutingAssayData" located in the
     * sprouting assay folder.
     * @throws IOException if errors
     */
    public void ExportData() throws IOException {
        Timestamp timestamp = new Timestamp(System.currentTimeMillis());
        ArrayList<Integer> percentages = new ArrayList<>();
        if (BATCH_RUN) {
            for (double percentage : HEPARIN_PERCENTAGES) {
                percentage = Math.round(percentage*100);
                percentages.add((int)percentage);
            }
        } else {
            percentages.add((int)Math.round(HEPARIN_PERCENTAGES[0]*100));
        }

        String timestamp_string = ((timestamp.toString().replace(" ","_").replace(".", "-").replace(":", "-")).substring(0, 10) +" " + (percentages) + "%");
        Path fileName= Path.of("SproutingAssay\\SproutingAssayData\\" + timestamp_string + "sensitivityall.csv");
        Path timeDataFileName =  Path.of("SproutingAssay\\SproutingAssayData\\" + timestamp_string + "sensitivityall_timeData.csv");
        Path headCellDistanceFileName = Path.of("SproutingAssay\\SproutingAssayData\\" + timestamp_string + "sensitivityall_headCellDistance.csv");
        int i = 1;
        while (Files.exists(fileName)){
            fileName= Path.of("SproutingAssay\\SproutingAssayData\\" + timestamp_string + " (" + i + ")" + "sensitivityallbaseline.csv");
            timeDataFileName =  Path.of("SproutingAssay\\SproutingAssayData\\" + timestamp_string + " (" + i + ")" +"sensitivityall_timeData.csv");
            headCellDistanceFileName = Path.of("SproutingAssay\\SproutingAssayData\\" + timestamp_string + " (" + i + ")" + "sensitivityall_headCellDistance.csv");
            i++;
        }
        StringBuilder dataset = CSV;
        Files.writeString(fileName, dataset);

        if (EXPORT_TIME_DATA){
            StringBuilder timeDataset = TIME_CSV;
            Files.writeString(timeDataFileName, timeDataset);
        }

        if (EXPORT_HEAD_CELL_DISTANCE_DATA){
            StringBuilder timeDataset = HEAD_CSV;
            Files.writeString(headCellDistanceFileName, timeDataset);
        }

    }

    /**
     * Initializes the headings for the CSV file.
     */
    public void Initialize_CSV(){
        CSV.append("Heparin Percentage (%), Total Vessel Length (microns), Fold Change (%), VEGF sensitivity, Diffusion Coefficient, VEGF Sensitivity, Initial percent head cells, VEGF intake, hours, vegfdegradation, low-med thres, med-high thres, vegf loaded, media exchange");
    }


    /**
     * Calculates and stores the fold change af intervals specified by FOLD_CHANGE_SAMPLE_TIME
     */
    public void CalculateFoldChangeOverTime(){
        if (GetTick()%FOLD_CHANGE_SAMPLE_TICKS == 0){
            FCTime.add(((double)GetTick())/TICKS_PER_HOUR);
            int currentVesselCount = countVessels();
            double foldChange = (double)(currentVesselCount)/initialCultureSize;
            foldChangeOverTime.add(foldChange);
        }
    }

    /**
     * Clears the foldChangeOverTime variable between model runs.
     */
    public void ClearFoldChange() {
        this.foldChangeOverTime.clear();
    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////      MAIN METHOD      /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    public static void main(String[] args) throws IOException {
        GridWindow gridWin = new GridWindow("Endothelial Cells",x, y, SCALE_FACTOR); // instantiates window for agents
        GridWindow VEGFWin = new GridWindow("VEGF Diffusion", x, y, SCALE_FACTOR); // instantiates window for diffusion

        sproutGrid model = new sproutGrid(x, y); // instantiates agent grid

        Path fileName= Path.of("SproutingAssay\\SproutingAssayData");  // sets Data export location
        File SproutingAssayDatafile = new File(String.valueOf(fileName));

        if (EXPORT_DATA && !SproutingAssayDatafile.exists()) {  // create a data export file if necessary
            if (!SproutingAssayDatafile.mkdir()) {
                throw new IOException("SproutingAssayData folder not made");
            }
        }

        if (!BATCH_RUN){ // if the run is not a batch run
            // Initialize
            model.initVesselsCircleCulture(model);  // initialize the spheroid
            model.initMAPParticles(model, HEPARIN_PERCENTAGES[0]); // initialize the MAP particles around the spheroid

            model.Initialize_CSV();  // initialize the CSV
            model.initialCultureSize= model.countVessels(); // keep track of the initial vessel count of spheroid (used in calculating fold change)
            for (int i = 0; i < TIMESTEPS; i++){ // while the model runs (for each tick)
                // pause
                gridWin.TickPause(TICK_PAUSE); // how fast the simulation runs
                // model step
                model.StepCells(); // prompt the cells to take action

                // draw
                model.DrawPDE(VEGFWin); // draw the PDE window
                model.DrawModel(gridWin); // draw the agent window
                model.CalculateFoldChangeOverTime(); // calculate fold change if necessary
            }
            if (EXPORT_DATA){ // if data export is set to TRUE,
                model.CollectData(HEPARIN_PERCENTAGES[0]); // then collect data
                model.ExportData(); // and export it.
            }

        } else {    // If this is a batch run,
            model.Initialize_CSV(); // initialie the CSV

            for (double heparinPercentage : HEPARIN_PERCENTAGES) { // For each percentage listed in HEPARIN_PERCENTAGES
                for (int trial = 0; trial < TRIALS; trial++) { // perform the amount of trials specified in TRIALS (for each trial...)
                    // initialize
                    model.Reset(); // reset the model for the next trial
                    model.ResetTick(); // reset the time tick for the next trial
                    model.ClearFoldChange(); // reset the fold change variable for the next trial
                    model.VEGF = new PDEGrid2D(x, y); // initialize the diffusion grid
                    model.initVesselsCircleCulture(model); // initialize spheroid
                    model.initMAPParticles(model, heparinPercentage); // initialize MAP particles

                    model.initialCultureSize= model.countVessels(); // keep track of the initial vessel count of spheroid (used in calculating fold change)
                    for (int i = 0; i < TIMESTEPS; i++){ // while the model runs (for each tick...)
                        // pause
                        gridWin.TickPause(TICK_PAUSE); // how fast the simulation runs
                        // model step
                        model.StepCells(); // prompt the cells to take action

                        // draw
                        model.DrawPDE(VEGFWin); // draw the PDE window
                        model.DrawModel(gridWin); // draw the agent window
                        model.CalculateFoldChangeOverTime(); // calculate fold change if necessary
                    }
                    if (EXPORT_DATA){ // if data export is set to TRUE,
                        model.CollectData(heparinPercentage); // then collect data
                    }
                }
            }
            if (EXPORT_DATA){  // If data export is TRUE,
                model.ExportData(); // Export all data when finished.
            }
        }
    }
}

