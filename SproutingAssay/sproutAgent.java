/*
GRIFFIN LAB
LAUREN PRUETT, ALEX TAING
UNIVERSITY OF VIRGINIA
2D Hybrid Model Sprouting Angiogenesis Publication 2022
*/

package SproutingAssay;

import HAL.GridsAndAgents.AgentSQ2D;
import HAL.Util;
import java.util.ArrayList;

public class sproutAgent extends AgentSQ2D<sproutGrid> {


    // DO NOT MODIFY FOR PARAMETERS
    public static int HEAD_CELL = 0;  // integer representing type HEAD cell
    public static int BODY_CELL = 1;  // integer representing type BODY cell
    public static int MAP_PARTICLE = 2;  // integer representing type MAP_PARTICLE
    public static int HEPARIN_MAP = 3;  // integer representing type HEPARIN_MAP particle

//    public static int HEAD_CELL_COLOR = Util.RED;
    public static int HEAD_CELL_COLOR = Util.GREEN; // head cells are green for visualization
    public static int BODY_CELL_COLOR = Util.RED; // body cells color
    public static int MAP_PARTICLE_COLOR = Util.RGB(128.0 / 255, 128.0 / 255, 128.0 / 255); // normal MAP;
    public static int HEPARIN_MAP_COLOR = Util.RGB(0.0 / 255, 0.0 / 255, 217.0 / 255); // Heparin MAP;

    // take all parameters for the sproutGrid class.
    public static double VESSEL_VEGF_INTAKE = sproutGrid.VESSEL_VEGF_INTAKE;
    public final static double VEGF_SENSITIVITY = sproutGrid.VEGF_SENSITIVITY;
    public final static int MAX_ELONGATION_LENGTH = sproutGrid.MAX_ELONGATION_LENGTH;
    public final static double PERSISTENCY_TIME = sproutGrid.PERSISTENCY_TIME;
    public final static int BRANCH_DELAY_TIME = sproutGrid.BRANCH_DELAY_TIME;
    public static final double HEP_MAP_VEGF_RELEASE = sproutGrid.HEP_MAP_VEGF_RELEASE;
    public static final double MEDIA_EXCHANGE_SCHEDULE_TICKS = sproutGrid.MEDIA_EXCHANGE_SCHEDULE_TICKS;
    public static final double VEGF_DEGRADATION_RATE = sproutGrid.VEGF_DEGRADATION_RATE;
    public final static double LOW_BRANCHING_PROBABILITY= sproutGrid.LOW_BRANCHING_PROBABILITY;
    public final static double LOW_MED_VEGF_THRESHOLD = sproutGrid.LOW_MED_VEGF_THRESHOLD;
    public final static double MED_BRANCHING_PROBABILITY= sproutGrid.MED_BRANCHING_PROBABILITY;
    public final static double MED_HIGH_VEGF_THRESHOLD = sproutGrid.MED_HIGH_VEGF_THRESHOLD;
    public final static double HIGH_BRANCHING_PROBABILITY= sproutGrid.HIGH_BRANCHING_PROBABILITY;
    public final static double REQUIRED_VEGF_GRADIENT_DIFFERENCE = sproutGrid.REQUIRED_VEGF_GRADIENT_DIFFERENCE;

    int color;  // color of cell in visualization
    int type;  // cell type (HEAD, BODY, MAP_PARTICLE, HEPARIN_MAP)
    int length = 0;  // length of the vessel (since starting location)
    int target; // simulates elongation direction:  cells will pursue a target coordinate in the direction of the highest VEGF concentration
    double migration_rate = sproutGrid.MIGRATION_RATE;  // endothelial cell migration rate
    double branching_probability; // probability that this cell will branch (set by BRANCHING_PROBABILITY in sproutGrid class)
    int since_last_growth;  // Works to realize a "migration rate".  A cell can only elongate one pixel every X ticks, determined by MIGRATION_RATE.  This stores the time since its last 1 pixel elongation.
    int elongation_length;  // threshold length which determines when an endothelial cell must find a new direction to migrate (based on highest VEGF concentration)
    int ticks_since_direction_change;  // // Works to realize a "persistency time".  A cell must persist in 1 direction for X time, determined by PERSISTENCY_TIME.  This stores the time since the cell's last redirection.
    int since_last_branch; // endothelial cells may branch to more frequently than BRANCH_DELAY_TIME.  This stores how long it has been since the vessel's last branch.


    /**
     * Gets the location with the highest VEGF concentration within the cell's radius of sight
     * @return returns the location of highest concentration of VEGF (if there are multiple highest, then it will return a random one of them)
     */
    public int HighestConcentrationVEGF() {
        assert G != null;
        int VEGF_options = G.MapEmptyHood(G.VEGFHood, Isq()); // gets the cell's range of detecting VEGF

        // Get the locations around the cell with highest concentration VEGF
        double maxConcentration = -1; // holds the max concentration so far (initially -1)
        ArrayList<Integer> maxConcentrationLocations = new ArrayList<>(); // holds the coordinates for the locations of highest concentration
        for (int i = 0; i < VEGF_options; i++) { // Iterates through all nearby coordinates
            double test_concentration = G.VEGF.Get(G.VEGFHood[i]); // gets the concentration at the nearby coordinate (called test_concentration)
            if ((test_concentration > maxConcentration) && (test_concentration > VEGF_SENSITIVITY)) { // if the concentration here is larger than the max so far (called maxConcentration),
                maxConcentration = test_concentration; // then set that concentration as the new max
                maxConcentrationLocations.clear(); // clear the old locations of highest concentration
                maxConcentrationLocations.add(G.VEGFHood[i]); // add this location to the list of highest concentrations
            } else if (test_concentration == maxConcentration) { // if the test_concentration is equal to the current max concentration
                maxConcentrationLocations.add(G.VEGFHood[i]); // add the coordinate to the running list of locations of highest concentration
            }
        }

        // if there were no highest concentrations, or if the concentration is 0, then return 0 (meaning highest VEGF coordinate was not found)
        if (maxConcentrationLocations.size() < 1) { // if there were no locations of highest concentration at all
            return 0;
        } else if (maxConcentration <= 0) { // if max concentration was 0
            return 0;
        }

        return maxConcentrationLocations.get((int) (Math.random() * maxConcentrationLocations.size())); // return a random one of the locations of highest concentration
    }

    /**
     * Given the location of a target as an int, this function will find best location for next cell duplication to move towards this location
     * (i.e. gives the next location for the cell to migrate towards VEGF, chemotaxis)
     *
     * @param target location in the direction of the highest VEGF concentration
     * @return int location that allows the vessel to grow closer to target (highest VEGF direction)
     */
    public int HoodClosestToTarget(int target) {

        int minDistance = Integer.MAX_VALUE; // holds the minimum distance of a coordinate to the target (gets updated with each location check)
        ArrayList<Integer> mincoordint = new ArrayList<>(); // keeps track of the candidates for the best migration location toward the target

        assert G != null;
        int options = G.MapEmptyHood(G.divHood, Isq()); // open areas around cell are found

        // find open areas (since cells can only migrate to empty areas)
        for (int i = 0; i < options; i++) { // iterate thorough the open areas
            int MAPcount = 0; // tally of how many MAP/HEP_MAP present
            for (sproutAgent cell : G.IterAgents(G.divHood[i])) { // iterate through all the cells at that coordinate
                if ((cell.type == MAP_PARTICLE) || (cell.type == HEPARIN_MAP)) { // if there is MAP/HEP_MAP there
                    MAPcount++; // then keep track that there was a particle there.
                    break; // and stop looking, since this space is occupied
                }
            }
            if (MAPcount == 0) { // If there were no occlusions with MAP particles, then check to see if it is close to the target point
                int[] hoodPoint = {G.ItoX(G.divHood[i]), G.ItoY(G.divHood[i])}; // this is the location that we are checking is/isn't a valid migration location

                // get the distance from the candidate migration point to target
                int dist = Math.abs((int) Math.hypot(G.ItoX(target) - hoodPoint[0], G.ItoY(target) - hoodPoint[1]));

                if (dist < minDistance) { // if the candidate point's distance is closer than the current minimum distance, then
                    minDistance = dist; // the minimum distance is updated to the new closest distance
                    mincoordint.clear();// the old list of candidates is cleared
                    mincoordint.add(G.I(hoodPoint[0], hoodPoint[1])); // and the new closest point is added to the empty list
                } else if (dist == minDistance) { // But, if the point is just as close as the ones on the list
                    mincoordint.add(G.I(hoodPoint[0], hoodPoint[1])); // it is added to the list of the candidate points that are just as close
                }
            }
        }
        if (mincoordint.size() == 0) { // if there are no best locations to migrate, then
            return 0; // return 0
        }
        return mincoordint.get((int) (Math.random() * mincoordint.size())); // otherwise, return a random point on the list of minimum distance coordinates
    }

    /**
     * Initializes a cell with color and type
     *
     * @param type type of cell/particle
     *
     */
    public void Init(int type) {

        this.type = type;

        if (type == HEAD_CELL) {
            this.color = HEAD_CELL_COLOR; // Growing vessel cells
        } else if (type == BODY_CELL) {
            this.color = BODY_CELL_COLOR; // vessel stalk cells
        } else if (type == MAP_PARTICLE) {
            this.color = MAP_PARTICLE_COLOR; // normal MAP
        } else if (type == HEPARIN_MAP) {
            this.color = HEPARIN_MAP_COLOR; // Heparin MAP
//            this.heparin_particle_release_VEGF = true;
        }
    }

    /**
     * Initializes a vessel cell.
     *
     * @param type type of cell/particle
     * @param length the current length of the vessel
     */
    public void InitVessel(int type, int length) {
        this.type = type;
        this.length = length;
        if (type == HEAD_CELL){
            this.color = HEAD_CELL_COLOR;
        } else if (type == BODY_CELL){
            this.color = BODY_CELL_COLOR;
        }
    }

    /**
     * Initializes a vessel cell with extra parameters allowing for implementation of migration rate
     *
     * @param type The celltype of the cell to be initialized
     * @param length The current length of the vessel up to this cell
     * @param target The target that the vessel is growing towards (highest concentration of VEGF)
     * @param since_last_elongation the time since the cell has grown a pixel (used to implement migration rate)
     * @param ticks_since_direction_change the amount of time since the vessel has last changed direction
     * @param elongation_length the current length since last direction change (current length traveled in current target direction)
     * @param since_last_branch the time since the last the the vessel has branched
     */
    public void InitVesselMigrationRate(int type, int length, int target, int since_last_elongation, int ticks_since_direction_change, int elongation_length, int since_last_branch) {
        this.type = type;
        this.length = length;
        this.target = target;
        this.since_last_growth = since_last_elongation;
        this.color = HEAD_CELL_COLOR;
        this.ticks_since_direction_change = ticks_since_direction_change;
        this.elongation_length = elongation_length;
        this.since_last_branch = since_last_branch;
    }

    /**
     * Calls vessel cells to consume appropriate amount of VEGF specified by VESSEL_VEGF_INTAKE
     */
    public void ConsumeVEGF() {
        assert G != null;
        if ((G.VEGF.Get(Isq()) >= 0) && ((type == HEAD_CELL) || (type == BODY_CELL))) { // Head cells and body cells consume VEGF
            G.VEGF.Add(Isq(), -VESSEL_VEGF_INTAKE); //Endothelial cell internalizes VEGF
            if (G.VEGF.Get(Isq()) < 0){ // if the VEGF concentration is less than 0,
                G.VEGF.Set(Isq(), 0); // then set it to zero.
            }
        }
    }

    /**
     * Replicates media exchange, simulating HEP_MAP sequestering and release of VEGF provided by fresh media, exchanged
     * as specificed by MEDIA_EXCHANGE_SCHEDULE
     */
    public void ExchangeMedia(){
        assert G != null;
        if (G.GetTick()%MEDIA_EXCHANGE_SCHEDULE_TICKS == 0){  // if the current tick is a multiple of  MEDIA_EXCHANGE_SCHEDULE_TICKS, then,
            if (type == HEPARIN_MAP){ // if the cell is of type HEPARIN_MAP
                G.VEGF.Set(Isq(), HEP_MAP_VEGF_RELEASE); // set the particle to release VEGF as prescribed by HEP_MAP_VEGF_RELEASE
            }
        }
    }


    /**
     * Performs vessel elongation analogous to that specified in Mehdizadeh et al.  The vessel will elongate towards the
     * highest concentration VEGF.  This function incorporates elongation, chemotaxis, and branching.  (see comments for
     * function details)
     */
    public void VesselGrowthByRate() {
        // checks of G is null (just for safely)
        assert G != null;

        // only continue for head cells
        if (type != HEAD_CELL){
            return;
        }

        // if the persistency time has not been reached, then increment it by 1
        if (ticks_since_direction_change <= PERSISTENCY_TIME) {
            ticks_since_direction_change += 1;
        }

        // if the time since last branch has not been reached, then increment it by 1
        if (since_last_branch <= BRANCH_DELAY_TIME){
            since_last_branch += 1;
        }

        // Makes the vessel grow at a certain rate (i.e. 30 microns/hr) (this is NOT elongation length!!!)
        if (since_last_growth < migration_rate){ // if it's not yet time to migrate, then don't migrate yet and stop any cell activity.
            since_last_growth += 1;
            return;
        }

        // only continue if there is enough VEGF around, if not then stay quiescent
        if (G.VEGF.Get(Isq()) < VEGF_SENSITIVITY) {
            return;
        }


        int cellDivLocation; // (stores the the location where the cell will grow to)

        // if persistency time has been reached, and the elongation length has been reached, then find a new direction to migrate.
        if ((( ticks_since_direction_change >= PERSISTENCY_TIME) && (elongation_length >= MAX_ELONGATION_LENGTH)) || (target == 0)){
            int highestConcentrationCoord = HighestConcentrationVEGF(); // find new target location: the location in the sight radius with the highest VEGF
            highestConcentrationCoord = CalculateTargetScaleTargetCoord(highestConcentrationCoord); // find a point in the same direction, but very far away, so you won't reach it: want to persist in that direction
            cellDivLocation = HoodClosestToTarget(highestConcentrationCoord); // and find the closest adjacent coordinate in the direction towards the highestConcentrationCoord
            if ((highestConcentrationCoord > 0) && (cellDivLocation > 0) && (G.PopAt(cellDivLocation) == 0)){ // If there is a location of highest concentration and there is an open adjacent spot... (the values will be 0 if none were found)
                if (!(G.VEGF.Get(Isq()) + REQUIRED_VEGF_GRADIENT_DIFFERENCE <= G.VEGF.Get(cellDivLocation))){ // check if there is enough difference in VEGF between current location and target growth spot, (if not, then end)
                    return;
                }
                G.NewAgentSQ(cellDivLocation).InitVesselMigrationRate(HEAD_CELL, this.length + 1, highestConcentrationCoord, since_last_growth, 0, 1, since_last_branch); // make a new head cell at cellDivLocation (vessel elongates)
                InitVessel(BODY_CELL, this.length); // and make the old cell a body cell
            }

            // branching
            if(G.rng.Double() < branching_probability && since_last_branch > BRANCH_DELAY_TIME){ // if the cell happens to branch (branching probability is a parameter of the cell, and is modified by a function CalculateBranchingProbability)
                // find the options for branching locations around the cell
                int options = MapEmptyHood(G.divHood);
                if (options >= 1) { // if there is an open nearby location, then branch there
                    cellDivLocation = G.divHood[G.rng.Int(options)]; // choose a random empty location
                    if (G.PopAt(cellDivLocation) == 0){  // if there are no agents at that location
                        if (!(G.VEGF.Get(Isq()) + REQUIRED_VEGF_GRADIENT_DIFFERENCE <= G.VEGF.Get(cellDivLocation))){ // if there is not a large enough difference of VEGF between the current location and the branching location, then do not branch. Otherwise,
                            return;
                        }
                        G.NewAgentSQ(cellDivLocation).InitVesselMigrationRate(HEAD_CELL, this.length + 1, 0, 0, this.ticks_since_direction_change =0, 1, 0); // make a new head cell at the branching location
                        InitVessel(BODY_CELL, this.length); // make the old cell a body cell
                    }
                }
            }
        } else if (elongation_length < MAX_ELONGATION_LENGTH){ // if the cell has not yet reached it's max elongation length, then continue to elongate
            cellDivLocation = HoodClosestToTarget(target); // find an open adjacent location closest to the target
            if (G.PopAt(cellDivLocation) == 0) {  // double check that there are no other agents at that location
                if (!(G.VEGF.Get(Isq()) + REQUIRED_VEGF_GRADIENT_DIFFERENCE <= G.VEGF.Get(cellDivLocation))){  // if the VEGF concentration difference between the current location and the elongation location is insufficient, then stop.  Else,
                    return;
                }
                G.NewAgentSQ(cellDivLocation).InitVesselMigrationRate(HEAD_CELL, this.length + 1, target, since_last_growth, ticks_since_direction_change, elongation_length + 1, since_last_branch); // make a new cell there
                InitVessel(BODY_CELL, this.length); // make the old cell a body cell
            }
        }
    }

    /**
     * This function takes in a target coordinate of a cell and views it in relation with that cell.  It then treats it
     * a vector, and scales the coordinate such that it is multiple times further, but still in the same direction
     * in relation to the cell.
     * @param targetCoord  the original coordinate target of a cell
     * @return a coordinate in the same direction, but scaled to be a larger distance than the original target.
     */
    public int CalculateTargetScaleTargetCoord(int targetCoord){
        assert G != null;
        int scale = 10; // initialize scale factor
        int target_x = G.ItoX(targetCoord); // the target X coordinate
        int target_y = G.ItoY(targetCoord); // the target y coordinate
        int current_x = Xsq(); // the current X coordinate of the cell
        int current_y = Ysq(); // the current y coordinate of the cell
        int scaled_vector_x = scale*(target_x-current_x); // takes the difference of the target and current x coordinates and scales the difference by the scale factor
        int scaled_vector_y = scale*(target_y-current_y); // takes the difference of the target and current y coordinates and scales the difference by the scale factor
        int[] new_target = {(scaled_vector_x+current_x),(scaled_vector_y+current_y)}; // reapplies the scaled to the current x and y coordinates
        while (!G.In(new_target[0], new_target[1])){ // if the scale is too large, then repeat the steps but for a smaller scale value
            scale -= 1;
            scaled_vector_x = scale*(target_x-current_x);
            scaled_vector_y = scale*(target_y-current_y);
            new_target[0] = scaled_vector_x+current_x;
            new_target[1] = (scaled_vector_y+current_y);
        }
        return G.I(new_target[0], new_target[1]); // return the scaled coordinate
    }

    /**
     *  Based on the current quantity of VEGF at the location of the cell, the relative branching probability is applied to the cell
     */
    public void CalculateBranchingProbability(){
        assert G != null;
        if (G.VEGF.Get(Isq()) < LOW_MED_VEGF_THRESHOLD){
            branching_probability = LOW_BRANCHING_PROBABILITY;
        } else if (G.VEGF.Get(Isq()) < MED_HIGH_VEGF_THRESHOLD){
            branching_probability = MED_BRANCHING_PROBABILITY;
        } else if (G.VEGF.Get(Isq()) > HIGH_BRANCHING_PROBABILITY){
            branching_probability = HIGH_BRANCHING_PROBABILITY;
        }
    }

    /**
     * Implements VEGF degradation such that at every 90 ticks, VEGF degrades at prescribed VEGF_DEGRADATION_RATE
     */
    public void VEGFDegrade() {
        assert G != null;
        if (G.VEGF.Get(Isq()) >= 0 && G.GetTick()%90 == 0)  { // VEGF degrades every 90 ticks consistent with half life of 90 min
            G.VEGF.Add(Isq(), -(G.VEGF.Get(Isq()))* VEGF_DEGRADATION_RATE);
            if (G.VEGF.Get(Isq()) < 0){
                G.VEGF.Set(Isq(), 0);
            }
        }
    }

    /**
     * Calls all functions for a single time tick.  This function is called in sproutGrid to initiate all cells to take
     * action at each tick.
     */
    public void StepCell() {
        
        // Eat VEGF
        ConsumeVEGF();

        // initialize VEGF
        ExchangeMedia();

        // calculate branching probability
        CalculateBranchingProbability();

        // Elongate
        VesselGrowthByRate();

        //Degrade
        VEGFDegrade();
    }
}
