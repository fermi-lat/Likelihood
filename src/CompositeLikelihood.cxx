//namespace Likleihood {

CompositeLikelihood::CompositeLikelihood() : optimizers::Statistic() {}

CompositeLikelihood::~CompositeLikelihood() throw() {}

void CompositeLikelihood::addComponent(const std::string & srcName, 
                                       LogLike & component) {
   m_components[srcName] = & component;
   
   m_normParIndex[srcName];
}

double CompositeLikelihood::value() const {
   double my_value(0);
   ComponentMap_t::ComponentConstIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
      my_value += it->second->value();
   }
   return my_value;
}

void CompositeLikelihood::
getFreeParams(std::vector<optimizers::Parameter> & params) const {
   if (m_components.empty()) {
      throw std::runtime_error("getFreeParams: empty composite list");
   }

   ComponentConstIterator_t it(m_components.begin());
   it->second->getFreeParams(params);
   ++it;
   for ( ; it != m_components.end(); ++it) {
      const Source & my_source(*it->second->getSource(it->first));
      params.push_back(
         const_cast<optimizers::Function &>(my_source.spectrum()).normPar());
   }
}

void CompositeLikelihood:setFreeParamValues(std::vector<double> & values) {
   size_t nsrcs(m_components.size());
//
// Set the parameter values for each component using the first
// component values, then reset the normalization parameters.
// 
   ComponentIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
      it->second->setFreeParamValues(values)
   }
   for (size_t i(nsrcs-1), it=m_components.end()-1; i > 0; --it, i--) {
      it->second->normPar().setValue(values.at(i));
   }
}

unsigned int CompositeLikelihood::getNumFreeParams() const {
   ComponentConstIterator_t it(m_components.begin());
   unsigned int npars(it->second->getNumFreeParams());
   return npars + m_components.size() - 1;
}

void CompositeLikelihood::getFreeDerivs(std::vector<double> & derivs) const {
   size_t nsrcs(m_components.size());
//
// Get the derivatives wrt to the free parameters of the first component.
// component values, then reset the normalization parameters.
// 
   ComponentConstIterator_t it(m_components.begin());
   for ( ; it != m_components.end(); ++it) {
      it->second->getFreeDerivs(derivs);
   }
   for (size_t i(nsrcs-1), it=m_components.end()-1; i > 0; --it, i--) {
      it->second->normPar().setValue(values.at(i));
   }
   
}

} // namespace Likleihood
