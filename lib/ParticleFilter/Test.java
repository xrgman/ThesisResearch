public class DoParticleFilter extends BasicFilter {
    public static final int amountOfParticles = 10000;
    private boolean particleFilterEnabled;

    public static double walkinglength = 75;
    private final List<Particle> allParticles = new ArrayList<>();
    private final List<Cell> allCells;
    private final List<Walls> allWalls;
    public Set<String> allowedCoordinates = new HashSet<>();

    public Set<String> allowedStairCoordinates = new HashSet<>();

    private int allowedPixels;
    private final int wallThickness;
    public static float stdevPosition = 8f;
    public static float stdevWalking = 10f;
    private final int noiseMEAN = 0;

    private final CustomLayoutView myCustomLayout;
    private Floor currentFloor;

    //Initialize all the lists and variables that are needed for the particle filter
    public DoParticleFilter(List<Cell> allCells, List<Walls> allWalls, int wallThickness, CustomLayoutView myCustomLayout){
        super(allCells);

        this.allCells = allCells;
        this.allWalls = allWalls;
        this.wallThickness = wallThickness;
        this.myCustomLayout = myCustomLayout;
        allowedPixels = 0;
        currentFloor = MainActivity.CurrentFloor;

        calculateAllowedCoordinates();
    }

    /**
     * Create particles according to bayes posterior
     * @param nrOfParticles - Nr. of particles to be plotted
     * @param posterior - The posterior to base the plotting of particles on.
     */
    public void createParticlesAccordingToPosterior(int nrOfParticles, List<Float> posterior) {
        //Clearing possible old particles:
        allParticles.clear();

        int particlesLeft = nrOfParticles;

        for(int i = 0; i < posterior.size(); i++) {
            Cell cell = allCells.get(i);
            int nrOfParticlesCell = Math.round(posterior.get(i) * nrOfParticles);

            //If there are not enough particles left:
            if(particlesLeft < nrOfParticlesCell) {
                nrOfParticlesCell = particlesLeft;
            }

            //Decreasing nr of particles that are left:
            particlesLeft -= nrOfParticlesCell;

            //Create all the particles:
            for(int y = 0; y < nrOfParticlesCell; y++){
                Set<String> coordinatesAllowed = cell.getAllowedCoordinates();

                allParticles.add(createParticle(coordinatesAllowed, cell.startxCoordinate + cell.width, cell.startYcoordinate + cell.height));
            }
        }

        myCustomLayout.updateParticles(allParticles, false);
    }
    
    public static float stdevWalking = 10f;
    private final int noiseMEAN = 0;

    //Add the new movement to the particles and add some noise onto it,
    //check if they are still valid else put them in the list of invalid once
    private void processNewMovement(double XMovement, double YMovement, double direction, boolean bayes){
        Random r = new Random(1);
        List<Particle> outOfBoundsParticles = new ArrayList<>();
        List<Particle> correctParticles = new ArrayList<>();
        int[] chooseCurrentCell = new  int[MainActivity.AMOUNTOFCELLS];    //used for choosing the current cell you are in
        Arrays.fill(chooseCurrentCell, 0);

        int noise1; int noise2; String coords;
        int sign1; int sign2;

        //Add some noise onto the values
        for(Particle particle : allParticles){
            do {
                sign1 = Math.random() >= 0.5 ?  -1 :  1;
                sign2 = Math.random() >= 0.5 ?  -1 :  1;
                noise1 = sign1 * (int)(r.nextGaussian() * stdevWalking + noiseMEAN);
                noise2 = sign2 * (int)(r.nextGaussian() * stdevWalking + noiseMEAN);
            } while (noise1 > walkinglength && noise2 > walkinglength);

            coords = (particle.xCoordinate + (int)XMovement + noise1) + "x" + (particle.yCoordinate + (int)YMovement + noise2);

            //Checking if coordinate is allowed:
            if(coordinateAllowedOnFloor(MainActivity.CurrentFloor, coords, bayes) &&
                    notThroughWall(particle.xCoordinate, particle.yCoordinate, particle.xCoordinate + (int)XMovement + noise1, particle.yCoordinate + (int)YMovement + noise2)) {
                particle.xCoordinate += (int)XMovement + noise1;
                particle.yCoordinate += (int)YMovement + noise2;
                correctParticles.add(particle);
                chooseCurrentCell[getCellOfParticle(particle)] += 1;
                //adjust weight of particle
            }else {
                outOfBoundsParticles.add(particle);
            }
        }

        Log.d("PARTICLE", "In total " + outOfBoundsParticles.size() + " particles were out of bound and " + correctParticles.size() + " were ok.");
        processOutOfBoundParticles(correctParticles, outOfBoundsParticles, chooseCurrentCell);
    }

    //WHAT TO DO IF ALL OUT OF BOUNDS??
    //order the list of correct particles based on their weight and then reset and add the outofboundsParticles to the
    //correctParticles
    private void processOutOfBoundParticles(List<Particle> correctParticles, List<Particle> outOfBoundsParticles, int[] chooseCurrentCell){
        EmpiricalDistribution fittedDistribution = new EmpiricalDistribution(correctParticles.size());  //bin count is correctParticles.size()
        double particleDistribution[] = new double[correctParticles.size()];


        //In the case that all the particles were out of bounds don't do anything and don't update
        if(outOfBoundsParticles.size() == amountOfParticles || correctParticles.size() == 0){
            //Check here for floor change, if floor has changed, spawn them inside that cell because that is probably where we are :)

            return;
        }

        //Increase the weight of all the correctparticles with 1

        for(int i = 0; i < correctParticles.size(); i++){
            correctParticles.get(i).setWeight(correctParticles.get(i).getWeight() + 1f/amountOfParticles);
            particleDistribution[i] = (double)correctParticles.get(i).getWeight();
        }

        fittedDistribution.load(particleDistribution);

        //correctParticles.sort(Comparator.comparing(Particle::weight));
        correctParticles.sort(Comparator.comparing(Particle::getWeight));
        Log.d("PARTICLE", "The highest particle x y coordinates are " + correctParticles.get(0).xCoordinate + " " + correctParticles.get(0).yCoordinate);
        double sumOfIncorrectParticles = 0;

        for(Particle wrongParticle : outOfBoundsParticles){
            //Choose randomly a new particle to fit onto, the higher weighted particles will be closer to the 0 thus faster chosen
            int index;
            Random r = new Random(1);
            do {
                //index = (int)(r.nextGaussian() * correctParticles.size()/5); //Need to choose this carefully
                index = (int)fittedDistribution.getNextValue();
            } while (index < 0 || index > correctParticles.size());

            //Update the likelihood of the chosen particle
            Particle chosenParticle = correctParticles.get(index);
            //chosenParticle.setWeight(chosenParticle.getWeight() + 1f/MainActivity.amountOfParticles);

            //Create some noise onto the x and y coordinates
            int noise1; int noise2; String coord;
            int newXCoord; int newYcoord;
            int sign1; int sign2;
            do {
                sign1 = Math.random() >= 0.5 ?  -1 :  1;
                sign2 = Math.random() >= 0.5 ?  -1 :  1;
                noise1 = sign1 * (int)(r.nextGaussian() * stdevPosition + noiseMEAN);
                noise2 = sign2 * (int)(r.nextGaussian() * stdevPosition + noiseMEAN);
                newXCoord = chosenParticle.xCoordinate + noise1;
                newYcoord = chosenParticle.yCoordinate + noise2;
                coord = newXCoord + "x" + newYcoord;
            } while (noise1 > walkinglength && noise2 > walkinglength && !allowedCoordinates.contains(coord));

            //Update the coordinates of the wrongParticle and reset the weight, NOtHING IS DONE YET WITH DIRECTION
            wrongParticle.setCoordinates(newXCoord, newYcoord, Floor.Second, chosenParticle.direction);
            wrongParticle.setWeight(1f/amountOfParticles);
            sumOfIncorrectParticles += 1f/amountOfParticles;
            chooseCurrentCell[getCellOfParticle(wrongParticle)] += 1;

        }

        normalizeWeightOfParticles(correctParticles, sumOfIncorrectParticles);

        //Just check which cell has the most amount of particles for coloring it
        myCustomLayout.setChosenCell(findCurrentCell(chooseCurrentCell));

        Log.d("PARTICLE", "Done with particle filter step");
    }

    private int findCurrentCell(int[] chooseCurrentCell) {
        int highestCount = 0;
        int cellNumber = 0;

        for(int i = 0; i < MainActivity.AMOUNTOFCELLS; i++){
            if(chooseCurrentCell[i] > highestCount){
                highestCount = chooseCurrentCell[i];
                cellNumber = i;
            }
        }

        Log.d("PARTICLE", "Highest cell is " + cellNumber + " with " + highestCount + " particles");

        float MINIMALPARTICLESCOUNT = amountOfParticles * (70f / 100f);

        if(highestCount > MINIMALPARTICLESCOUNT){
            return cellNumber;
        }

        return -1;
    }

    //Add all the weight together, the ones from incorrect particles are given as parameter and then normalize
    private void normalizeWeightOfParticles(List<Particle> correctParticles, double sumOfIncorrectParticles){

        for (Particle particle : correctParticles){
            sumOfIncorrectParticles += particle.getWeight();
        }

        for(Particle particle : correctParticles){
            particle.setWeight((float)(particle.getWeight()/sumOfIncorrectParticles));
        }
    }

    //Method that compares two lines (the wall and the particle movement line) and checks if they intersect
    private boolean notThroughWall(int startX, int startY, int finishX, int finishY){
        Line2D line1 = new Line2D( startX,  startY,  finishX,  finishY);
        Line2D line2;
        for(Walls wall : allWalls){
            //Skip walls not on current floor:
            if(wall.floor != MainActivity.CurrentFloor) {
                continue;
            }

            line2 = new Line2D(wall.startXCoordinate, wall.startYCoordinate, wall.stopxCoordinate, wall.stopYCoordinate);
            if(Line2D.intersects(line1, line2)){
                return false;
            }
        }
        return true;
    }
}