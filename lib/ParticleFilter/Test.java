
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


    //WHAT TO DO IF ALL OUT OF BOUNDS??
    //order the list of correct particles based on their weight and then reset and add the outofboundsParticles to the
    //correctParticles
 

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
