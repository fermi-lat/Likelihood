#ifndef Likelihood_UnbinnedLikelihood_h
#define Likelihood_UnbinnedLikelihood_h

namespace Likelihood {

class UnbinnedLikleihood : public LikelihoodBase {

public:

   UnbinnedLikelihood(const Observation & observation);

   ~UnbinnedLikelihood();

   double value() const;

   void getFreeDerivs(std::vector<double> & derivs) const;

   /// Used by UnbinnedAnalysis
   void computeEventResponses();
   void set_ebounds(double emin, double emax);

};

} // namespace Likelihood

#endif // Likelihood_UnbinnedLikelihood_h
