/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
 *             David Strelak (davidstrelak@gmail.com)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "reconstruction/movie_alignment_deformation_model.h"

void ProgMovieAlignmentDeformationModel::readParams()
{
	fnMovie = getParam("-i");
	fnMicrograph = getParam("-o");
    initDose = getDoubleParam("--initDose");
    perFrameDose = getDoubleParam("--perFrameDose");
    maxIterations = getIntParam("--maxIterations");
    upScaling = getIntParam("--upscaling");
    fnUnaligned = getParam("--ounaligned");
    fnGlobAligned = getParam("--oGlobAligned");
    fnGlobFrames = getParam("--oGlobFrames");
    fnOutFrames = getParam("--oCorrFrames");
    threadNumbers = getIntParam("-j");
    fnDark = getParam("--dark");
    fnGain = getParam("--gain");
    shiftLimit = getDoubleParam("--shiftLimit");
    patchShiftLimit = getDoubleParam("--patchShiftLimit");
    filterCutOff = getDoubleParam("--filterCut");
    downsample = getDoubleParam("--downsample");
    if (patchShiftLimit < 0) {
        patchShiftLimit = shiftLimit;
    }
    show();
}

void ProgMovieAlignmentDeformationModel::show()
{
    if (!verbose)
        return;
    std::cout 
    << "Input movie:          " << fnMovie           << std::endl
    << "Output micrograph:    " << fnMicrograph      << std::endl
    << "Initial dose:         " << initDose          << std::endl
    << "Per frame dose:       " << perFrameDose      << std::endl
    << "Max iterations:       " << maxIterations     << std::endl
    << "Up scaling coef.:     " << upScaling         << std::endl
	<< "Unaligned micrograph: " << fnUnaligned       << std::endl
    << "Globally aligned:     " << fnGlobAligned     << std::endl
    << "Globally corr. frames:" << fnGlobFrames      << std::endl
    << "Corrected frames:     " << fnOutFrames       << std::endl
    << "Threads number:       " << threadNumbers     << std::endl
    << "Dark image:           " << fnDark            << std::endl
    << "Gain image:           " << fnGain            << std::endl
    << "Shift limit:          " << shiftLimit        << std::endl
    << "Patch shift limit:    " << patchShiftLimit   << std::endl
    << "Filter cut off freq:  " << filterCutOff      << std::endl
    << "Downsamplin scale:    " << downsample        << std::endl;
}

void ProgMovieAlignmentDeformationModel::defineParams()
{
    addUsageLine("Align a set of frames by cross-correlation of the frames");
    addParamsLine("   -i <metadata>               : Metadata with the list of frames to align");
    addParamsLine("   -o <fn=\"\"> 		          : Give the name of a micrograph to generate an aligned micrograph");
    addParamsLine("  [--initDose <s=0>]           : Radiation dose received before first frame is taken");
    addParamsLine("  [--perFrameDose <s=0>]       : Radiation dose received after imaging each frame");
    addParamsLine("  [--maxIterations <N=20>]	  : Number of robust least squares iterations");
    addParamsLine("  [--upscaling <s=1>]          : UpScaling coefficient for super resolution image generated from model application");
    addParamsLine("  [--ounaligned <fn=\"\">]     : Give the name of a micrograph to generate an unaligned (initial) micrograph");
    addParamsLine("  [--oGlobAligned <fn=\"\">]   : Path to micrograph calculated from frames after just global correction");
    addParamsLine("  [--oGlobFrames <fn=\"\">]    : Where to save individual frames after global alignment");
    addParamsLine("  [--oCorrFrames <fn=\"\">]    : Where to save individual motion corrected frames");
    addParamsLine("  [-j <N=5>]                   : Maximum threads the program is allowed to use");
    addParamsLine("  [--dark <fn=\"\">]           : Dark correction image");
    addParamsLine("  [--gain <fn=\"\">]           : Gain correction image");
    addParamsLine("  [--shiftLimit <s=200>]       : Limits the maximal shift global alignment can calculate");
    addParamsLine("  [--patchShiftLimit <s=-1>]    : Limits the maximal shift alignment of patches can calculate (if -1 'shiftLimit' value is used");
    addParamsLine("  [--filterCut <s=0.05>]       : Low pass filter cut off requency");
    addParamsLine("  [--downsample <s=2.0>]       : Downsampling scale");
}

void ProgMovieAlignmentDeformationModel::run()
{  
    loadMovie(fnMovie, frames, timeStamps, fnDark, fnGain);
    int imageCount = frames.size();

    // --Initialize data structures
    globalShiftsX.resize(imageCount);
    globalShiftsY.resize(imageCount);

    partitions.resize(PARTITION_COUNT * PARTITION_COUNT);
    for (int i = 0; i < partitions.size(); i++) {
        partitions[i].resize(imageCount);
    }

    localShiftsX.resize(imageCount * PARTITION_COUNT * PARTITION_COUNT);
    localShiftsY.resize(imageCount * PARTITION_COUNT * PARTITION_COUNT);

    deformationCoeffsX.setlength(9);
    deformationCoeffsY.setlength(9);

    double MAX_SHIFT_THRESHOLD = 0.1;

    // --Calculations
    if (!fnUnaligned.isEmpty()) {  // save unaligned average (if selected)
        MultidimArray<double> unalignedMicrograph;
        averageFrames(frames, unalignedMicrograph);
        saveMicrograph(fnUnaligned, unalignedMicrograph);
    }

    std::cout << "Estimating global shifts" << std::endl;
    estimateShifts(frames, globalShiftsX, globalShiftsY, maxIterations,
            MAX_SHIFT_THRESHOLD, shiftLimit);
    std::cout << "Found global shifts (x, y):" << std::endl;
    for (int i = 0; i < globalShiftsX.size(); i++) {
        std::cout << i << ": (" << globalShiftsX[i] << ", " 
            << globalShiftsY[i] << ")" << std::endl;
    }
    std::cout << "Applying global shifts" << std::endl;
    applyShifts(frames, globalShiftsX, globalShiftsY);

    if (!fnGlobFrames.isEmpty()) { // save individual globally corrected frames
        FileName fn;
        for (size_t i = 0; i < frames.size(); i++) {
            fn.compose(fnGlobFrames + "Glob", i + 1, "mrc");
            std::cout << fn << std::endl;
            saveMicrograph(fn, frames[i]);
        }
    }

    if (!fnGlobAligned.isEmpty()) { // save globaly algined average
        MultidimArray<double> tmp;
        averageFrames(frames, tmp);
        saveMicrograph(fnGlobAligned, tmp);
    }
    std::cout << "Partitioning" << std::endl;
    partitionFrames(frames, partitions, PARTITION_COUNT);
    std::cout << "Estimating local shifts" << std::endl;
    estimateLocalShifts(partitions, localShiftsX, localShiftsY, maxIterations,
            MAX_SHIFT_THRESHOLD, patchShiftLimit);
    
    std::cout << "Estimating deformation model coefficients" << std::endl;
    calculateModelCoefficients(localShiftsX, timeStamps,
            deformationCoeffsX, frames[0].ydim, frames[0].xdim);
    calculateModelCoefficients(localShiftsY, timeStamps,
            deformationCoeffsY, frames[0].ydim, frames[0].xdim);

    std::cout << "Applying local motion correction" << std::endl;
    motionCorrect(frames, timeStamps, deformationCoeffsX, deformationCoeffsY,
            upScaling);

    if (!fnOutFrames.isEmpty()) { // save individual corrected frames
        FileName fn;
        for (size_t i = 0; i < frames.size(); i++) {
            fn.compose(fnOutFrames, i + 1, "mrc");
            std::cout << fn << std::endl;
            saveMicrograph(fn, frames[i]);
        }
    }

    std::cout << "Saving resulting average" << std::endl;
    MultidimArray<double> result;
    averageFrames(frames, result);
    saveMicrograph(fnMicrograph, result);
}

void ProgMovieAlignmentDeformationModel::calcLPF(
        double targetOccupancy, const MultidimArray<double>& lpf) {
    double iNewXdim = 1.0 / 4096;
    double sigma = targetOccupancy / 6; 
    double K = -0.5 / (sigma * sigma);
    for (int i = STARTINGX(lpf); i <= FINISHINGX(lpf); ++i) {
        double w = i * iNewXdim;
        A1D_ELEM(lpf, i) = exp(K * (w * w));
    }
}

void ProgMovieAlignmentDeformationModel::loadMovie(FileName fnMovie,
        std::vector<MultidimArray<double> >& frames,
        std::vector<double>& timeStamps, FileName fnDark, FileName fnGain)
{
    Image<double> dark, gain;
    if (not fnDark.isEmpty()) {
        dark.read(fnDark);
    }
    if (not fnGain.isEmpty()) {
        gain.read(fnGain);
        gain() = 1.0 / gain();
        double avg = gain().computeAvg();
        if (std::isinf(avg) || std::isnan(avg)) {
            REPORT_ERROR(ERR_ARG_INCORRECT,
                "The input gain image is incorrect, its inverse produces infinite or nan");
        }
    }

    Image<double> movieStack;
    movieStack.read(fnMovie);
    size_t Xdim, Ydim, Zdim, Ndim;
    movieStack.getDimensions(Xdim, Ydim, Zdim, Ndim);
    std::cout << "Loading stack with dimensions (x, y, z, n): ";
    std::cout << Xdim << "," << Ydim << "," << Zdim <<"," << Ndim << std::endl;
    bool switched = false;
    if (Zdim == 1 and Ndim > 1) {
        Zdim = Ndim;
        switched = true;
    }
    frames.resize(Zdim, MultidimArray<double>(Ydim, Xdim));
    for (size_t z = 0; z < Zdim; z++) {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(frames[z]) {
            if (not switched) {
                DIRECT_A2D_ELEM(frames[z], i, j) = DIRECT_NZYX_ELEM(movieStack(),
                                                                0, z, i, j);
            } else {
                DIRECT_A2D_ELEM(frames[z], i, j) = DIRECT_NZYX_ELEM(movieStack(),
                                                                z, 0, i, j);
            }
        }

        if (dark().xdim > 0) {
            frames[z] -= dark();
        }
        if (gain().xdim > 0) {
            frames[z] *= gain();
        }
        frames[z].setXmippOrigin();

        timeStamps.push_back(initDose + perFrameDose * z);
    }
}

void ProgMovieAlignmentDeformationModel::scaleLPF(
        const MultidimArray<double>& lpf, int xSize, int ySize,
        double targetOccupancy, MultidimArray<double>& result) {
    Matrix1D<double> w(2);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(result)
    {
        FFT_IDX2DIGFREQ(i, ySize, YY(w));
        FFT_IDX2DIGFREQ(j, xSize, XX(w));
        double wabs = w.module();
        if (wabs <= targetOccupancy)
            A2D_ELEM(result, i, j) = lpf.interpolatedElement1D(wabs * xSize);
    }
}

void ProgMovieAlignmentDeformationModel::filterAndBinFrame(
        MultidimArray<double>& frame, FourierFilter& filter, int xdim, int ydim)
{
    selfScaleToSizeFourier(xdim, ydim, frame, threadNumbers);
    filter.applyMaskSpace(frame);
}

void ProgMovieAlignmentDeformationModel::estimateShifts(
        std::vector<MultidimArray<double> >& data,
		std::vector<double>& shiftsX, std::vector<double>& shiftsY,
        int maxIterations, double minShiftTermination, double maxShift)
{
    int xdim = data[0].xdim;
    int ydim = data[0].ydim;
    std::vector<MultidimArray<double> > filtData;
    filtData.resize(data.size(), MultidimArray<double>(ydim, xdim));
    for (int t = 0; t < filtData.size(); t++) {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(data[t]) {
            DIRECT_A2D_ELEM(filtData[t], i, j) = DIRECT_NZYX_ELEM(data[t],
                                                            0, 0, i, j);
        }
    }

    FourierFilter filter;
    filter.FilterBand=LOWPASS;
    filter.FilterShape=RAISED_COSINE;
    filter.w1=filterCutOff;
    filter.generateMask(filtData[0]);

    for (int i = 0; i < filtData.size(); i++) {
        filtData[i].setXmippOrigin();
        filterAndBinFrame(filtData[i], filter, xdim / 2, ydim / 2);
    }
    std::cout << "Filter applied" << std::endl;

    std::vector<MultidimArray<double> > shiftedData;
    shiftedData.resize(data.size());
    for (int i = 0; i < shiftedData.size(); i++) {
        shiftedData[i].setXmippOrigin();
        filtData[i].setXmippOrigin();
        shiftedData[i].initZeros(filtData[0]);
        shiftedData[i] += filtData[i];
    }

	// prepare sum of images
    MultidimArray<double> sum;
    averageFrames(shiftedData, sum);

    saveMicrograph("/scratch/workdir/1.jpg", filtData[0]);
    saveMicrograph("/scratch/workdir/2.jpg", filtData[1]);
    saveMicrograph("/scratch/workdir/3.jpg", filtData[2]);

    // estimate the shifts
    double shiftX, shiftY;
    CorrelationAux aux;
    for (int cycle = 0; cycle < maxIterations; cycle++) {
        std::cout << "Cycle: " << cycle << std::endl;
        double topDiff = 0;
        for (int i = 0; i < data.size(); i++) {
            sum -= shiftedData[i];
            bestShift(sum, filtData[i], shiftX, shiftY, aux, NULL, maxShift);
            
            topDiff = std::max(std::max(std::abs(shiftsY[i] - shiftY),
                        std::abs(shiftsX[i] - shiftX)),
                    topDiff);

            shiftsY[i] = shiftY;
            shiftsX[i] = shiftX;

            //if (i > 0) {
                //if (i == 1) {
                    //translate(BSPLINE3, shiftedData[0], filtData[0],
                            //vectorR2(shiftsX[i], shiftsY[i]));
                //}
                translate(BSPLINE3, shiftedData[i], filtData[i],
                        vectorR2(shiftsX[i], shiftsY[i]));
                sum += shiftedData[i];
            //}
        }

        // recalculate shifts so, first frame has shift 0
        //for (int i = shiftsX.size() - 1; i >= 0; i--) {
            //shiftsX[i] -= shiftsX[0];
            //shiftsY[i] -= shiftsY[0];
        //}

        std::cout << "- top diff: " << topDiff << std::endl;

        std::vector<double> kX = {0, 0, -5, -5, -5, 2, -12, -10, -10};
        std::vector<double> kY = {0, -10, 2, -12, -5, -5, -5, 0, -10};
        double sumDiff = 0;
        for (int i = 0; i < shiftsX.size(); i++) {
            sumDiff += std::abs(kX[i] + (shiftsX[i] - shiftsX[0]));
            sumDiff += std::abs(kY[i] + (shiftsY[i] - shiftsY[0]));
        }
        std::cout << "- diff metic: " << sumDiff << ";  " << std::endl;

        std::cout << "- shifts: ";
        for (int i = 0; i < shiftsX.size(); i++) {
            double sx = shiftsX[i];
            double sy = shiftsY[i];
            std::cout << "(" << sx << "," << sy << "), ";
        }
        std::cout << std::endl;

        if (topDiff <= minShiftTermination) {
            // if further calculated shifts are only minimal end calculation
            std::cout << "Shift threshold reached: " << topDiff << std::endl;
            break;
        }

        // recalculate avrg
        for (int i = 0; i < data.size(); i++) {
        }
        averageFrames(shiftedData, sum);
    }

    // recalculate shifts so, first frame has shift 0
    for (int i = shiftsX.size() - 1; i >= 0; i--) {
        shiftsX[i] -= shiftsX[0];
        shiftsX[i] *= 2;
        shiftsY[i] -= shiftsY[0];
        shiftsY[i] *= 2;
    }
}

void ProgMovieAlignmentDeformationModel::estimateLocalShifts(
        std::vector<std::vector<MultidimArray<double> > >& partitions,
        std::vector<double>& shiftsX, std::vector<double>& shiftsY,
        int maxIterations, double minShiftTermination, double maxShift)
{
    //shiftsX and shiftsY contains shifts for all partitions []
    //shifts are organized 
	int partsPerFrame = partitions.size();
	int partDepth = partitions[0].size(); //frame count
    std::vector<double> tmpXShifts(partDepth);
    std::vector<double> tmpYShifts(partDepth);
    for (int i = 0; i < partitions.size(); i++) {
        std::cout << "Local movement estimation for partition " << i
            << std::endl;
        estimateShifts(partitions[i], tmpXShifts, tmpYShifts, maxIterations,
                minShiftTermination, maxShift);
        for (int j = 0; j < partDepth; j++) {
        	shiftsX[i + j*partsPerFrame] = tmpXShifts[j];
        	shiftsY[i + j*partsPerFrame] = tmpYShifts[j];
        }
    }
}

void ProgMovieAlignmentDeformationModel::applyShifts(
        std::vector<MultidimArray<double> >& data,
        const std::vector<double>& shiftsX, const std::vector<double>& shiftsY)
{
	for (int i = 0; i < data.size(); i++) {
		selfTranslate(BSPLINE3, data[i], vectorR2(shiftsX[i], shiftsY[i]),
                false, 0.0);
	}
}

void ProgMovieAlignmentDeformationModel::partitionFrames(
        const std::vector<MultidimArray<double> >& frames,
        std::vector<std::vector<MultidimArray<double> > >& partitions,
        int edgeCount)
{
    // properly resize the individual partitions
    int partSizeY = YSIZE(frames[0]) / edgeCount;
    int partSizeX = XSIZE(frames[0]) / edgeCount;
    int yReminder = YSIZE(frames[0]) - (partSizeY * edgeCount);
    int xReminder = XSIZE(frames[0]) - (partSizeX * edgeCount);

    for (int i = 0; i < partitions.size(); i++) {
        int partY = i / edgeCount;
        int partX = i % edgeCount;
        int xSize = partSizeX + (partX < xReminder ? 1 : 0);
        int ySize = partSizeY + (partY < yReminder ? 1 : 0);
        for (int j = 0; j < partitions[i].size(); j++) {
            partitions[i][j].resize(1, 1, ySize, xSize);
            partitions[i][j].initZeros();
            partitions[i][j].setXmippOrigin();
        }
    }

    int longerPartY = yReminder * (partSizeY + 1);
    int longerPartX = xReminder * (partSizeX + 1);
    for (int fi = 0; fi < frames.size(); fi++) {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(frames[fi]) {
            int partY, partX; // index of the partition
            int innerY, innerX; // index inside the partition
            
            if (i < longerPartY) { // partitions are longer along y axis here
                partY = i / (partSizeY + 1);
                innerY = i % (partSizeY + 1);
            } else {
                partY = yReminder + (i - longerPartY) / partSizeY;
                innerY = (i - longerPartY) % partSizeY;
            }

            if (j < longerPartX) { // partitions are longer along x axis here
                partX = j / (partSizeX + 1);
                innerX = j % (partSizeX + 1);
            } else {
                partX = xReminder + (j - longerPartX) / partSizeX;
                innerX = (j - longerPartX) % partSizeX;
            }
            dAij(partitions[partY * edgeCount + partX][fi], innerY, innerX) = dAij(frames[fi], i, j);
        }
    }
}

void ProgMovieAlignmentDeformationModel::calculatePartitionSize(int partIndex,
        int edgeCount, int frameHeight, int frameWidth, int& partXSize,
        int& partYSize)
{
	partYSize = frameHeight / edgeCount;
    partXSize = frameWidth / edgeCount;
    int yReminder = frameHeight - (partYSize * edgeCount);
    int xReminder = frameWidth - (partXSize * edgeCount);

    int partIndexY = partIndex / edgeCount;
    int partIndexX = partIndex % edgeCount;
    partYSize = partYSize + (partIndexY < yReminder ? 1 : 0);
    partXSize = partXSize + (partIndexX < xReminder ? 1 : 0);
}

void ProgMovieAlignmentDeformationModel::calculateModelCoefficients(
        const std::vector<double>& shifts,
        const std::vector<double>& timeStamps, alglib::real_1d_array& c,
        int frameHeight, int frameWidth)
{
	for (size_t i = 0; i < c.length(); i++) {
		c[i] = 0.05;  // define initial guess
	}

	alglib::real_1d_array y;
	y.setlength(shifts.size());
	for (size_t i = 0; i < y.length(); i++) {
		y[i] = shifts[i];
	}

	alglib::real_2d_array positions;
	positions.setlength(shifts.size(), 3);
	int partsInFrame = PARTITION_COUNT * PARTITION_COUNT;
    double cummulativeX, cummulativeY;
    int partSizeX, partSizeY;
	for (size_t i = 0; i < positions.rows(); i++) {
		int frameIndex = i / partsInFrame;
		int partIndex = i % partsInFrame;
        if (partIndex == 0) { //starting next frame
            cummulativeY = 0;
            cummulativeX = 0;
        } else if (partIndex % PARTITION_COUNT == 0) { //new partition line
            cummulativeX = 0;
            cummulativeY += partSizeY;
        }
        calculatePartitionSize(partIndex, PARTITION_COUNT, frameHeight,
                frameWidth, partSizeX, partSizeY);

		positions[i][0] = cummulativeY + (partSizeY / 2.0);
		positions[i][1] = cummulativeX + (partSizeX / 2.0);
		positions[i][2] = timeStamps[frameIndex];

        cummulativeX += partSizeX;
	}

    alglib::ae_int_t info;
    alglib::lsfitstate state;
    alglib::lsfitreport rep;
    double epsf = 0.000001;
    double epsx = 0.000001;
    alglib::ae_int_t maxits = 0;
    double diffstep = 0.000001;

    alglib::lsfitcreatef(positions, y, c, diffstep, state);
    alglib::lsfitsetcond(state, epsf, epsx, maxits);
    alglib::lsfitfit(state, pixelShiftAlg);
    alglib::lsfitresults(state, info, c, rep);
    //TODO: check state and saner default values
}

void ProgMovieAlignmentDeformationModel::pixelShiftAlg(
        const alglib::real_1d_array &c, const alglib::real_1d_array &dim,
        double &func, void *ptr)
{
	double y = dim[0];
	double x = dim[1];
	double t = dim[2];
    func = pixelShift(x, y, t, c);
}

double ProgMovieAlignmentDeformationModel::pixelShift(double x, double y, 
        double t, const alglib::real_1d_array& c)
{
    return (c[0] + c[1]*x + c[2]*x*x + c[3]*y + c[4]*y*y + c[5]*x*y) * 
    		(c[6]*t + c[7]*t*t + c[8]*t*t*t);
}

void ProgMovieAlignmentDeformationModel::motionCorrect(
        std::vector<MultidimArray<double> >& data,
        const std::vector<double>& timeStamps, const alglib::real_1d_array& cx,
        const alglib::real_1d_array& cy, int scaling)
{
    int origWidth = data[0].xdim;
    int origHeight = data[0].ydim;
    MultidimArray<double> tmp(origHeight * scaling, origWidth * scaling);
	for (int i = 0; i < data.size(); i++) {
        std::cout << "..." << i << std::endl;
        revertDeformation(data[i], tmp, cx, cy, timeStamps[i], scaling);
        scaleToSize(BSPLINE3, data[i], tmp, origWidth, origHeight);
    }
}

void ProgMovieAlignmentDeformationModel::revertDeformation(
        MultidimArray<double>& input, MultidimArray<double>& output,
        const alglib::real_1d_array& cx, const alglib::real_1d_array& cy,
		double t, int scalingFactor)
{   
    input.resetOrigin();
    output.resetOrigin();
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(output)
    {
        int y = i;
        int x = j;

        int nonScaledY = y / scalingFactor;
        int nonScaledX = x / scalingFactor;

        double posy = nonScaledY - pixelShift(nonScaledX, nonScaledY, t, cy);
        double posx = nonScaledX - pixelShift(nonScaledX, nonScaledY, t, cx);

        double val = input.interpolatedElement2DOutsideZero(posx, posy);
        dAij(output, y, x) = val;
    }
    //input.resetOrigin();
    //output.resetOrigin();
}

void ProgMovieAlignmentDeformationModel::averageFrames(
        const std::vector<MultidimArray<double> >& data,
        MultidimArray<double>& out)
{
	out.initZeros(data[0]);
    out.setXmippOrigin();
    for (int i = 0; i < data.size(); i++) {
        out += data[i];
    }
}

void ProgMovieAlignmentDeformationModel::saveMicrograph(
        const FileName fnMicrograph, const MultidimArray<double>& micrograph)
{
	Image<double> img(micrograph);
	img.write(fnMicrograph);
}

