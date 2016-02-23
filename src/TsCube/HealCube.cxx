/**
 * @file TsCube.cxx
 * @brief Application for producing "test-statistic" cubes.
 *
 * @author E. Charles, from TsMap by J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/src/TsCube/TsCube.cxx,v 1.3 2015/12/10 00:59:24 echarles Exp $
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <sstream>

#include "facilities/commonUtilities.h"

#include "astro/HealpixProj.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "st_facilities/Environment.h"
#include "st_facilities/Util.h"

#include "evtbin/Binner.h"

#include "optimizers/dArg.h"
#include "optimizers/Optimizer.h"
#include "optimizers/OptimizerFactory.h"
#include "optimizers/Exception.h"

#include "Likelihood/AppHelpers.h"
#include "Likelihood/BinnedLikelihood.h"
#include "Likelihood/LogLike.h"
#include "Likelihood/SourceMap.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/ScanUtils.h"
#include "Likelihood/FitScanner.h"

using namespace Likelihood;

/**
 * @class TsMap
 * @brief Class for encapsulating methods for creating a test-statistic
 * map.
 *
 */

class HealCube : public st_app::StApp {
public:
  HealCube();
  virtual ~HealCube() throw() {
    try {
      delete m_opt;
      delete m_helper;
    } catch (std::exception & eObj) {
      std::cerr << eObj.what() << std::endl;
    } catch (...) {
    }
  }
  virtual void run();
  virtual void banner() const;
private:

  AppHelpers * m_helper;
  st_app::AppParGroup & m_pars;
  Likelihood::BinnedLikelihood * m_logLike;
  optimizers::Optimizer* m_opt;
  st_stream::StreamFormatter * m_formatter;
  Likelihood::FitScanner* m_scanner;
  astro::HealpixProj* m_proj;
  astro::SkyDir m_refDir; 
  std::string m_coordSys;

  void promptForParameters();
  void readSrcModel();
  void selectOptimizer();
  void setGrid();
  void computeMap();
  void writeFitsFile();
  void setPointSourceSpectrum(PointSource &src);
  
  static std::string s_cvs_id;
};

st_app::StAppFactory<HealCube> myAppFactory("gthealcube");

HealCube::HealCube() 
  : st_app::StApp(), m_helper(0), 
    m_pars(st_app::StApp::getParGroup("gthealcube")),
    m_logLike(0), m_opt(0),
    m_formatter(new st_stream::StreamFormatter("gthealcube", "", 2)),
    m_scanner(0),
    m_proj(0),
    m_refDir(){
  setVersion(s_cvs_id);
}

std::string HealCube::s_cvs_id("$Name:  $");

void HealCube::banner() const {
  int verbosity = m_pars["chatter"];
  if (verbosity > 2) {
    st_app::StApp::banner();
  }
}

void HealCube::promptForParameters() {
  m_pars.Prompt("cmap");
  m_pars.Prompt("bexpmap");   
  m_pars.Prompt("expcube");
  m_pars.Prompt("irfs");
  m_pars.Prompt("srcmdl");
  m_pars.Prompt("outfile");
  m_pars.Prompt("hpx_ordering_scheme");
  m_pars.Prompt("hpx_order");
  m_pars.Prompt("hpx_region");
  m_pars.Prompt("coordsys");
  m_pars.Prompt("xref");
  m_pars.Prompt("yref");

  m_pars.Save();
}

void HealCube::run() {
  promptForParameters();
  
  m_helper = new AppHelpers(&m_pars, "BINNED");
  
  m_helper->checkOutputFile();
  std::string expcube = m_pars["expcube"];
  std::string irfs = m_pars["irfs"];
  if (expcube != "" && expcube != "none") {
    m_helper->observation().expCube().readExposureCube(expcube);
  }
  
  std::string cmap = m_pars["cmap"];
  m_helper->setRoi(cmap, "", false);
  if (!m_helper->observation().expCube().haveFile()) {
    throw std::runtime_error
      ("An exposure cube file is required for binned analysis. "
       "Please specify an exposure cube file.");
  }
  st_facilities::Util::file_ok(cmap);
  CountsMap* dataMap = new CountsMap(cmap);
  bool apply_psf_corrections = m_pars["psfcorr"];
  bool computePointSources(true);
  m_logLike = new BinnedLikelihood(*dataMap, 
				   m_helper->observation(),
				   cmap, 
				   computePointSources,
				   apply_psf_corrections);
  std::string bexpmap = m_pars["bexpmap"];
  AppHelpers::checkExposureMap(cmap, bexpmap);
  dynamic_cast<BinnedLikelihood *>(m_logLike)->setVerbose(false);
  
  readSrcModel();
  selectOptimizer();
  setGrid();
  computeMap();
  writeFitsFile();
}


void HealCube::readSrcModel() {
  std::string srcModelFile = m_pars["srcmdl"];
  if (srcModelFile != "" && srcModelFile != "none") {
    st_facilities::Util::file_ok(srcModelFile);
    bool requireExposure = false;
    bool loadMaps = false;
    bool addPointSources = true;
    m_logLike->readXml(srcModelFile, m_helper->funcFactory(),
		       requireExposure, addPointSources,
		       loadMaps);
  }
}

void HealCube::selectOptimizer() {
  std::string optimizer = m_pars["optimizer"];
  m_opt = optimizers::OptimizerFactory::instance().create(optimizer,
							  *m_logLike);
}

void HealCube::setGrid() {

  int hpx_order = m_pars["hpx_order"];
  std::string hpx_scheme = m_pars["hpx_ordering_scheme"];
  std::string hpx_region = m_pars["hpx_region"];

  double xref = m_pars["xref"];
  double yref = m_pars["yref"];
  std::string coordsys = m_pars["coordsys"];
  bool is_galactic(coordsys == "GAL");
  int nnorm = m_pars["nnorm"];
  std::string testSourceName = m_pars["target"];

  m_refDir() = astro::SkyDir(xref,yref,
			     (is_galactic ? astro::SkyDir::GALACTIC : astro::SkyDir::EQUATORIAL ) )();   


  Healpix_Ordering_Scheme scheme = hpx_scheme == "RING" ? RING : NEST;

  astro::HealpixProj* hpxProj = new astro::HealpixProj(hpx_order,scheme,is_galactic);
  m_scanner = new FitScanner(*m_logLike,*m_opt,*hpxProj,hpx_region);
  m_proj = hpxProj;

  int status = testSourceName.empty() ? 
    m_scanner->setPowerlawPointTestSource(m_helper->funcFactory()) : 
    m_scanner->setTestSourceByName(testSourceName);
  if ( status ) {
    throw std::runtime_error("Failed to make a powerlaw point test source");
    return;
  }  
}


void HealCube::computeMap() {

  // Get optimizer options
  int verbosity = m_pars["chatter"];
  verbosity -= 2;
  double tol = m_pars["ftol"];
  std::string tol_type = m_pars["toltype"];
  optimizers::TOLTYPE tolType(optimizers::ABSOLUTE);
  if (tol_type == "REL") {
      tolType = optimizers::RELATIVE;
  }

  // These is always true for this application
  static const bool doTsMap(true);
  static const bool doSED(true);

  // Hidden parameters for the loop
  int nnorm = m_pars["nnorm"];
  int maxiter = m_pars["maxiter"];
  int ST_scan_level = m_pars["stlevel"];
  bool remakeTestSource = m_pars["remakesrc"];
  double normSigma = m_pars["nsigma"];
  double covScale_bb = m_pars["covscale_bb"]; 
  double covScale = m_pars["covscale"]; 

  //m_scanner->set_verbose_null(3);
  //m_scanner->set_verbose_bb(3);
  //m_scanner->set_verbose_scan(3);

  int status = m_scanner->run_tscube(doTsMap,doSED,nnorm,normSigma,covScale_bb,covScale,
				     tol,maxiter,tolType,remakeTestSource,ST_scan_level);
  
  if ( status != 0 ) {
    // Go ahead and let it try to write the output, for now
    // throw std::runtime_error("FitScanner::run_tscube() returned with error");
    return;
  }
  return;
}


void HealCube::writeFitsFile() {
  std::string data_dir = facilities::commonUtilities::getDataPath("Likelihood");
  std::string template_file = facilities::commonUtilities::joinPath(data_dir, "TsCubeTemplate");
  int status = m_scanner->writeFitsFile(m_pars["outfile"],
					"gthealcube",
					template_file);
  return;
}
