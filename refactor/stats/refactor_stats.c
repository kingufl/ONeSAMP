// Most of this code is verbatim from the original OneSamp
#include "refactor_stats.h"


/*! \def fallingQuotient(s, t1, t2, c)
 *  \brief Returns s * ((t1)/(t2)) * ((t1-1)/(t2-1)) * ... ((t1-c+1)/(t2-c+1))
 *  The value c must be nonnegative.
 */
double fallingQuotient(double s, double t1, double t2, int c){
  return (c == 0) ? s : fallingQuotient(s * t1 / t2, t1 - 1, t2 - 1, c - 1);
}

/*! \def allelePr()
 *  \brief Returns likelihood of finding val1 alleles in one place and val2 in another
 */
double allelePr(int val1, int val2, double theta){
  int n = val1 + val2;
  double result = fallingQuotient(1, n, theta + n - 1, n);
  int j;
  if(val1 != 0) result *= theta / (double) val1;
  if(val2 != 0) result *= theta / (double) val2;
  if(val1 == val2) result /= 2;
  return result;
}

/*! \def minAlleleCount()
 *  \brief Returns minimum number of mutated alleles in coalescent sample
 */
int minAlleleCount(){
  return 4;
}

/*! \def mutate(ALLELE_TYPE *gene)
 *  \brief Mutates the given gene based on whether it is a microsat or SNP
 */
void mutate(ALLELE_TYPE *gene, int motif){
  if(parseFormFlag() != 1) { mutateSNP(gene); } else { mutateMicroSat(gene, motif); }
}

/*! \def mutateSNP(ALLELE_TYPE *gene)
 *  \brief Mutates a SNP
 */
void mutateSNP(ALLELE_TYPE *gene){
  *gene = 2 * randBit() + randBit() + 1;
}

/*! \def mutateMicroSat(ALLELE_TYPE *gene)
 *  \brief Mutates a microsatellite
 */
void mutateMicroSat(ALLELE_TYPE *gene, int motif){
  if(*gene != 0) *gene += randBit() ? motif : -motif; // Random increase or decrease in microsat length.
  if(*gene < 12) *gene = 12;
  if(*gene > 996) *gene = 996;
}

/*! \def assort(int next_gen_count, struct gtype_type *offvec, struct gtype_type *mothers, struct gtype_type *fathers, int current_gen_count)
 *  \brief Generates the next generation of individuals based on the current one
 */
void assort(int next_gen_count, struct gtype_type *offvec, struct gtype_type *mothers, struct gtype_type *fathers, int current_gen_count, int samp, int num_loci)
{
  // Read in mutation rate
  int pMutation = parseMRate(samp);
  // Read in number of loci
  int j,i,m,d;
  // Simulate each individual's genotype
  for(j = 0; j < next_gen_count; ++j){
    // Select random mother from current generation
    m = disrand(0,current_gen_count-1);
    // Select random father from current generation
    d = disrand(0,current_gen_count-1);
    // For each allele, select one from parent and possibly mutate
    for(i = 0; i < num_loci; ++i){
      // Generate a new individual
      offvec[j].mgtype[i] = randBit() ? mothers[m].mgtype[i] : mothers[m].pgtype[i];
      if(branchWithProbability(pMutation)) mutate(offvec[j].mgtype + i, parseFormFlag() != 1 ? 0 : getMotifLengths()[i]);
      // Generate a new individual
      offvec[j].pgtype[i] = randBit() ? fathers[d].mgtype[i] : fathers[d].pgtype[i];
      if(branchWithProbability(pMutation)) mutate(offvec[j].pgtype + i, parseFormFlag() != 1 ? 0 : getMotifLengths()[i]);
    }
  }
}

/*! \def sortM(int **numberOfAlleles, int num_samples, double m[], int ***gType, int ****gcountPtr)
 *  \brief computes allele length divided by allele length range (doesn't work for SNPs)
 */
void sortM(int **numberOfAlleles, int num_samples, double m[], int ***gType, int ****gcountPtr)
{
  int ***gcount = *gcountPtr;
  int numloci = parseNLoci();
  int i,h,iloc,r,s,skip,mono;
  double Msum,M;
  int samp;
  for(samp=0;samp<num_samples;++samp){
    if(parseFormFlag() != 1){m[samp] = -1; continue; }
    Msum = 0.0;
    M = 0.0;
    mono = 0;
    for(iloc=0;iloc<numloci;++iloc){
      skip = 0;
      r = 1;
      s = 0;
        for(i=0;i<(numberOfAlleles[samp][iloc]-1);++i){
          for(h=i+1;h<numberOfAlleles[samp][iloc];++h){
            s = abs(gType[samp][iloc][i] - gType[samp][iloc][h]) + 1;
            if ((r < s) && (gcount[samp][iloc][i] > 0) && (gcount[samp][iloc][h] > 0)) r = s;
          }
        }
        for(i=0;i<numberOfAlleles[samp][iloc];++i) {
          if(gcount[samp][iloc][i] == 0) ++skip;
        }
      if(r==1) {++mono;continue;}
      Msum += (double)(numberOfAlleles[samp][iloc] - skip) / r;
    }
    if(mono == numloci) {m[samp] = 0.0; continue;}
    m[samp] = Msum /((double) numloci - mono);
  }
}

// STATISTIC 2: iis: Burrows/Weir estimator

/*! \def twolocusiis(int **numberOfAlleles, int num_samples, struct gtype_type **samp_data, double iis[], int ***gType, int ****gcountPtr)
 *  \brief computes composite LD estimator with alleles
 */
void twolocusiisAssist(int **numberOfAlleles, int num_samples, struct gtype_type **samp_data, double iis[], int ***gType, int ****gcountPtr, int samp)
{
 // if(0 == 0) iis[samp] = -1;
 // if(0 == 0) return;
  int ***gcount = *gcountPtr;
  // Total number of allele pairs in this sample
  int prs = 0; 
  // Accumulate this statistic for each sample
  iis[samp] = 0;
  // Index of locus i
  int iloc;
  // Keep track of progress through the loops (for debug output)
  int progress = 0;

  int jloc;
  int alprs;
  int al1;
  int al2;
  double psq1obs;
  double psq2obs;
  int cnt;
  int ofreq1;
  int ofreq2;
  int ind;
  int doublesum;
  gtype_type curIndex;
  int miData;
  int mjData;
  int piData;
  int pjData;
  int gi1;
  int gj2;
  double dcnt;
  double p1;
  double p2;
  double jointAB;
  double d1;
  double d2;
  double sqrtFactor1;
  double sqrtFactor2;
  double r_intermediate;
  double r;
  double sampCorrection;
  double result = 0;
  #pragma omp parallel for private(jloc, alprs, al1, al2, psq1obs, psq2obs, cnt, ofreq1, ofreq2, ind, doublesum, curIndex, miData, mjData, piData, pjData, gi1, gj2, dcnt, p1, p2, jointAB, d1, d2, sqrtFactor1, sqrtFactor2, r_intermediate, r, sampCorrection)
  for(iloc = 0; iloc < parseNLoci(); iloc++){
    // Index of locus j
    int jloc;
    for(jloc = iloc + 1; jloc < parseNLoci(); jloc++){
      progress++;
      // Count of number of allele pairs
      int alprs = 0;
      // Index of allele 1
      int al1;
      for(al1 = 0; al1 < numberOfAlleles[samp][iloc]; al1++){
        // Index of allele 2
        int al2;
        for(al2 = 0; al2 < numberOfAlleles[samp][jloc]; al2++){
          // Observed count of homozygotes in first locus
          psq1obs = 0;
          // Observed count of homozygotes in second locus
          psq2obs = 0;
          // Count of frequency of indices with both entries specified
          cnt = 0;
          // Observed frequency of allele 1
          ofreq1 = 0;
          // Observed frequency of allele 2
          ofreq2 = 0;
          // Index into current sample
          // ind;
          // Count of matches pairs of genotypes
          doublesum = 0;
          for(ind = 0; ind < parseInputSamples(); ind++){
            curIndex = samp_data[samp][ind];
            miData = curIndex.mgtype[iloc];
            mjData = curIndex.mgtype[jloc];
            piData = curIndex.pgtype[iloc];
            pjData = curIndex.pgtype[jloc];
            gi1 = gType[samp][iloc][al1];
            gj2 = gType[samp][jloc][al2];
            #define miMatch (miData == gi1)
            #define mjMatch (mjData == gj2)
            #define piMatch (piData == gi1)
            #define pjMatch (pjData == gj2)
            #define iANDmatch (miMatch && piMatch)
            #define jANDmatch (mjMatch && pjMatch)
            #define iORmatch  (miMatch || piMatch)
            #define jORmatch  (mjMatch || pjMatch)
            if(miData != 0 && mjData != 0){
              cnt++;
              if(miData == gi1 && piData == gi1) {++psq1obs; ofreq1 += 2;}
              else if ((miData == gi1 && piData != gi1) || (miData != gi1 && piData == gi1)) {++ofreq1;}
              if      (mjData == gj2 && pjData == gj2) {++psq2obs; ofreq2 += 2;}
              else if ((mjData == gj2 && pjData != gj2) || (mjData != gj2 && pjData == gj2)) {++ofreq2;}
              if (iANDmatch && jANDmatch) {doublesum += 4; continue;}
              if (iANDmatch && jORmatch)  {doublesum += 2; continue;}
              if (iORmatch && jANDmatch)  {doublesum += 2; continue;}
              if (iORmatch && jORmatch)   {doublesum++; continue;}
            }
          }
          if(ofreq1 == 2 * cnt) continue; // Uncomment this to skip monoallelic sites
          if(ofreq2 == 2 * cnt) continue; // Uncomment this to skip monoallelic sites
          // Double value of cnt
          dcnt = (double) cnt;
          // p1 frequency in computation
          p1 = ofreq1 / (2 * dcnt);
          // p2 frequency in computation
          p2 = ofreq2 / (2 * dcnt);
          // Joint frequency
          jointAB = doublesum / (2 * dcnt);
          // Departure from Hardy-Weinberg equilibrium
          d1 = psq1obs / dcnt - p1 * p1;
          // Departure from Hardy-Weinberg equilibrium
          d2 = psq2obs / dcnt - p2 * p2;
          // Factor in denominator of r
          sqrtFactor1 = sqrt(p1 * (1 - p1) + d1);
          // Factor in denominator of r
          sqrtFactor2 = sqrt(p2 * (1 - p2) + d2);
          // Value of r
          // R intermediate
          r_intermediate = jointAB - 2 * p1 * p2;
          // Value of r
          r = r_intermediate / sqrtFactor1 / sqrtFactor2;
          // Correction for having a finite sample
          sampCorrection = dcnt / (dcnt - 1);
          // Final value of r
          r *= sampCorrection;
          // Accumulate the sum of the r squared values

          // If r is not NaN
          if(!(r == r)){
            // Disequilibrium in data
            // Backtrack and follow LDNe's strategy for computing if we get NaN
            d1 = 0;
            d2 = 0;
            sqrtFactor1 = sqrt(p1 * (1 - p1) + d1);
            // Factor in denominator of r
            sqrtFactor2 = sqrt(p2 * (1 - p2) + d2);
            // Value of r
            // R intermediate
            r_intermediate = jointAB - 2 * p1 * p2;
            // Value of r
            r = r_intermediate / sqrtFactor1 / sqrtFactor2;
            // Correction for having a finite sample
            sampCorrection = dcnt / (dcnt - 1);
            // Final value of r
            r *= sampCorrection;
            continue;
          } 
          // Debug output
          //printf("cnt = %d\n", cnt);
          //printf("psq1obs = %f\npsq2obs = %f\n", psq1obs, psq2obs);
          //printf("ofreq1 = %f\nofreq2 = %f\n", ofreq1, ofreq2);
          //printf("al1 = %d\nal2 = %d\n", al1, al2);
          //printf("iloc = %d\njloc = %d\np1 = %f\np2 = %f\n", iloc, jloc, p1, p2);
          //printf("d1 = %f\nd2 = %f\n", d1, d2);
          //printf("sqrtFactor1 = %f\nsqrtFactor2 = %f\n", sqrtFactor1, sqrtFactor2);
          //printf("jointAB = %f\n", jointAB);
          //printf("r_intermediate = %f\n", r_intermediate);
          //printf("r = %f\n", r);
          //printf("iis[samp] (current) = %f\n\n", iis[samp]);
          #pragma omp atomic
          result += r*r;
          #pragma omp atomic
          alprs++;
        }
      }

      // Debug output
      // int priorPercentage = 200 * (progress - 1) / (parseNLoci() * (parseNLoci() - 1));
      // int currentPercentage = 200 * progress / (parseNLoci() * (parseNLoci() - 1));
      // if(currentPercentage != priorPercentage) printf("Progress computing statistic iis (sample no. %d): %d%\n", samp + 1, currentPercentage);
      #pragma omp atomic
      prs += alprs;
    }
  }
  // Take the average of the r squared values
  result /= prs;
  iis[samp] = result;
}

void twolocusiis(int **numberOfAlleles, int num_samples, struct gtype_type **samp_data, double iis[], int ***gType, int ****gcountPtr)
{
  // Index of sample
  int samp;
  // Iterate through each sample
  for(samp = 0; samp < num_samples; samp++) {
    twolocusiisAssist(numberOfAlleles, num_samples, samp_data, iis, gType, gcountPtr, samp);
  }
}

// STATISTIC 3: lnbeta: imbalance in allele lengths

/*! \def beta(int **numberOfAlleles, int num_samples, double lnbeta[], int ***gType, int ****gcountPtr)
 *  \brief computes imbalance in allele lengths (no SNPs)
 */
void beta(int **numberOfAlleles, int num_samples, double lnbeta[], int ***gType, int ****gcountPtr)
{
  int final_indivs_count = parseInputSamples();
  int ***gcount = *gcountPtr;
  int numloci = parseNLoci();
  int samp,iloc,kal,skip;
  double psq,sumlen,meanlen,po,varlen,cvarlen,cvarpo,beta;
  for(samp=0;samp<num_samples;++samp) {
    if(parseFormFlag() != 1) {lnbeta[samp] = -1; continue;}
    beta = 0.0;
    skip = 0;
    for(iloc=0;iloc<numloci;++iloc) {
      sumlen = psq = cvarpo = cvarlen = varlen = 0.0;

      for(kal=0;kal<numberOfAlleles[samp][iloc];++kal) {
        if (gcount[samp][iloc][kal] == 0) continue;
        psq += ((double)gcount[samp][iloc][kal]*gcount[samp][iloc][kal])/((double)4 * final_indivs_count * final_indivs_count);
        sumlen += gcount[samp][iloc][kal]*gType[samp][iloc][kal];
      }

      if (psq == 1.0) {++skip; continue;}

      meanlen = (sumlen/(2 * final_indivs_count));
      po = (psq * 2 * final_indivs_count-1)/(2 * final_indivs_count-1);

      for(kal=0;kal<numberOfAlleles[samp][iloc];++kal) {
      varlen += gcount[samp][iloc][kal]*((gType[samp][iloc][kal]-meanlen)*(gType[samp][iloc][kal]-meanlen));
      }

      cvarlen = 2.0*varlen/(2*final_indivs_count-1);
      cvarpo = ((1/(po*po))-1)/2.0;
      beta += (log(cvarlen)-log(cvarpo));
    }
  if(numloci == skip) {lnbeta[samp] = 0.0; continue;}
  lnbeta[samp] = beta/(numloci-skip);
  }
}

// STATISTIC 4: hetx: Wright's Fis: SNPs okay, formula is unchanged
// STATISTIC 5: mnehet: expected mean heterozygosity: SNPs okay, formula is unchanged
// The unbiased sample statistic is from the Nei 1987 source.

/*! \def hetexcess(int **numberOfAlleles,int num_samples, struct gtype_type **samp_data, double *hetx, double *mnehet, int ***gType, int ****gcountPtr)
 *  \brief computes excess heterozygosity statistics
 */
void hetexcess(int **numberOfAlleles,int num_samples, struct gtype_type **samp_data, double *hetx, double *mnehet, int ***gType, int ****gcountPtr)
{
  int ***gcount = *gcountPtr;
  int numloci = parseNLoci();
  int final_indivs_count = parseInputSamples();
  int iloc, al1, dblp, ind, skipind, nonzeroindices, skiploc;
  double sumhobs, sumhexp, obshomo, exphomo;

  int samp;
  for(samp = 0; samp < num_samples; samp++){
    sumhobs = sumhexp = skiploc = 0;
    for(iloc = 0; iloc < numloci; iloc++) {
      obshomo = exphomo = nonzeroindices = 0;
      dblp = 0;
      skipind = 0;
      for(ind = 0; ind < final_indivs_count; ind++) {
        if((samp_data[samp][ind].mgtype[iloc] == 0) || (samp_data[samp][ind].pgtype[iloc] == 0)) ++skipind;
        else if(samp_data[samp][ind].mgtype[iloc] == samp_data[samp][ind].pgtype[iloc]) ++dblp; // Count of homozygotes
      }
      ind -= skipind; // ind now counts number of legal pairs
      obshomo += (double)dblp / ind; // Frequency of homozygotes
      for(al1 = 0; al1 < numberOfAlleles[samp][iloc]; al1++) {
        if(gType[samp][iloc][al1] == 0){ continue; }
        nonzeroindices += gcount[samp][iloc][al1];
        exphomo += gcount[samp][iloc][al1] * gcount[samp][iloc][al1]; // Expected frequency of homozygotes
      } // als per locus

      if(nonzeroindices == 0 || ind == 0 || numberOfAlleles[samp][iloc] == 1 || (numberOfAlleles[samp][iloc] == 2 && (gType[samp][iloc][0] == 0 || gType[samp][iloc][1] == 0))) { skiploc++; continue; } // Number of heterozygotes is undefined or monoallelic site.
      exphomo /= nonzeroindices * nonzeroindices;
      double observedHeterozygoteFrequency = 1 - obshomo;
      double expectedHeterozygoteFrequency = 1 - exphomo;
      // printf("iloc = %d\nobservedH = %f\n", iloc, observedHeterozygoteFrequency);
      // printf("expectedH = %f\n",expectedHeterozygoteFrequency);
      double sampleCorrectionFactor = ((double) ind) / (ind - 1);
      sumhobs += observedHeterozygoteFrequency; // Accumulate actual frequencies of heterozygotes

      double samplehexp = sampleCorrectionFactor * (expectedHeterozygoteFrequency - observedHeterozygoteFrequency/(2*ind));

      // Bounds check
      if(!(samplehexp > 0)) samplehexp = 0;
      if(!(samplehexp < 1)) samplehexp = 1;

      sumhexp += samplehexp; // Accumulate expected frequencies of heterozygotes
    } // loci
    mnehet[samp] = sumhexp / (double) (numloci - skiploc);

    //if(mnehet[samp] != mnehet[samp]) {
    //  printf("Illegal mean NEHET\n"); writeoutput(samp_data);
    //}

    hetx[samp] = (sumhexp == 0) ? 1/0.0 : 1 - sumhobs / sumhexp;

    //if(hetx[samp] != hetx[samp]) {
    //  printf("Illegal excess het\n"); writeoutput(samp_data);
    //}

  } // samples
}

// STATISTIC 6: mnals: Compute this statistic first.

/*! \def counts(int ***numberOfAllelesPtr, struct gtype_type **samp_data, double mnals[], int ***gType, int ****gcountPtr)
 *  \brief Generates genotype counts and mean number of allele data
 */
void counts(int ***numberOfAllelesPtr, struct gtype_type **samp_data, double mnals[], int ***gType, int ****gcountPtr)
{
  int locusID; // Identity of locus
  int i; // Loop index
  int **numberOfAlleles = *numberOfAllelesPtr;
  int ***gcount = *gcountPtr;

  int samp;
  // For each sample
  for(samp = 0; samp < parseIterations(); samp++){
    int p = 0;

    // For each locus
    for(locusID = 0; locusID < parseNLoci(); ++locusID){
      numberOfAlleles[samp][locusID] = 0;

      // Iterate through each individual genes through each pair
      // Apply a counting sort algorithm
      int j;
      for(j = 0; j < 2 * parseInputSamples();++j) {
        // Current individual
        int indiv = j / 2;
        // Current value of allele genotype
        int val = (j % 2 == 1) ? samp_data[samp][indiv].pgtype[locusID] : samp_data[samp][indiv].mgtype[locusID];
        for(i = 0; i < numberOfAlleles[samp][locusID];++i){
          if(val == gType[samp][locusID][i]){
            ++gcount[samp][locusID][i];
            break;
          }
        }

        // Create a new index to store information for a new allele
        if(i == numberOfAlleles[samp][locusID]){
          gcount[samp][locusID][i] = 1;
          gType[samp][locusID][i] = val;
          ++numberOfAlleles[samp][locusID];
        }
      }

      // Do this only if we have more than one allele at this locus
      // ???
      if(numberOfAlleles[samp][locusID] != 1){
        for(i = 0; i < parseInputSamples(); ++i){
          samp_data[samp][i].pgtype[p] = samp_data[samp][i].pgtype[locusID];
          samp_data[samp][i].mgtype[p] = samp_data[samp][i].mgtype[locusID];
        }
        numberOfAlleles[samp][p] = numberOfAlleles[samp][locusID];
        for(i = 0; i < numberOfAlleles[samp][p]; i++) {
          gType[samp][p][i] = gType[samp][locusID][i];
          gcount[samp][p][i] = gcount[samp][locusID][i];
        }
        p++;
      }
    }

    // Accumulate counts in numberOfAlleles data structure and store in mnals
    mnals[samp] = 0;
    for(locusID = 0; locusID < parseNLoci(); locusID++){
      mnals[samp] += numberOfAlleles[samp][locusID];
    }
    mnals[samp] /= parseNLoci();
  }
}

// FIXME unsure of licensing for this code? This was the original source for it.

// calculate first four moments of multilocus homozygosity, taken from p613 Numerical Receipes in C

// STATISTIC 7: mhomo: mean of multi-locus homozygosity: SNPs okay, formula is unchanged
// STATISTIC 8: varhomo: variance of multi-locus homozygosity : SNPs okay, formula is unchanged
// STATISTIC 9: skew
// STATISTIC 10: kurtosis

/*! \def multih(int num_samples, struct gtype_type **samp_data, double mhomo[], double varhomo[], double skhomo[], double kurhomo[], int ***gType)
 *  \brief computes moments of homozygosity
 */
void multih(int num_samples, struct gtype_type **samp_data, double mhomo[], double varhomo[], double skhomo[], double kurhomo[], int ***gType)
{
  int final_indivs_count = parseInputSamples();
  int samp, ind, i, cnt;
  int *data = (int *) malloc(final_indivs_count * sizeof(int));
  double s, ep, p, sdev;
  for(samp = 0; samp < num_samples; samp++)  {
    s = 0;
    for(ind = 0; ind < final_indivs_count; ind++) {
      cnt = 0;
      for(i = 0; i < parseNLoci(); i++)  {
        if(samp_data[samp][ind].mgtype[i] == samp_data[samp][ind].pgtype[i])  ++cnt;
      }
      data[ind] = cnt;
      s += cnt;
    }

    // printf("%f%d\n", s, final_indivs_count);
    mhomo[samp] = s/(double)final_indivs_count;

    sdev = ep = varhomo[samp] = skhomo[samp] = kurhomo[samp] = 0.0;
    for(i = 0; i < final_indivs_count; i++) {
      s = data[i] - mhomo[samp];
      ep += s;
      varhomo[samp] += (p = s*s);
      skhomo[samp] += (p *= s);
      kurhomo[samp] += (p *= s);
    }

    varhomo[samp] = (varhomo[samp]-ep*ep/final_indivs_count)/(final_indivs_count-1);
    sdev = sqrt(varhomo[samp]);
    if (varhomo[samp]) {
      skhomo[samp] /= (final_indivs_count*varhomo[samp]*sdev);
      kurhomo[samp] /= (final_indivs_count*varhomo[samp]*varhomo[samp]);
      kurhomo[samp] -= 3.0;
    }
  } // samp

  free(data);
}
