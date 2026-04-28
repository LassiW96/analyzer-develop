///////////////////////////////////////////////////////////////////////////////
//
//   THaScalerEvtHandler
//
//   Event handler for Hall A scalers.
//   R. Michaels,  Sept, 2014
//
//   This class does the following
//      For a particular set of event types (here, event type 140)
//      decode the scalers and put some variables into global variables.
//      The global variables can then appear in the Podd output tree T.
//      In addition, a tree "TS" is created by this class; it contains
//      just the scaler data by itself.  Note, the "fName" is concatenated
//      with "TS" to ensure the tree is unique; further, "fName" is
//      concatenated with the name of the global variables, for uniqueness.
//      The list of global variables and how they are tied to the
//      scaler module and channels is in the scaler.map file, or could
//      be hardcoded here.
//      NOTE: if you don't have the scaler map file (e.g. db_LeftScalevt.dat)
//      there will be no variable output to the Trees.
//
//   To use in the analyzer, your setup script needs something like this
//       gHaEvtHandlers->Add (new THaScalerEvtHandler("Left","HA scaler event type 140"));
//
///////////////////////////////////////////////////////////////////////////////

#include "ScalerEvtHandler.h"
#include "CodaDecoder.h"        // for CodaDecoder, THaEvData, coda_format_error
#include "Database.h"           // for CFile, VarType, ChopPrefix
#include "Module.h"             // for Module
#include "TClass.h"             // for TClass
#include "TError.h"             // for Error
#include "THaAnalysisObject.h"  // for THaAnalysisObject
#include "THaCrateMap.h"        // for THaCrateMap
#include "THaGlobals.h"         // for gHaVars
#include "THaVar.h"             // for THaVar
#include "THaVarList.h"         // for THaVarList
#include "TROOT.h"              // for TROOT
#include "TString.h"            // for TString, operator==, operator+
#include "TTree.h"              // for TTree
#include "Textvars.h"           // for Trim, vsplit
#include <algorithm>            // for find_if, fill_n
#include <cassert>              // for assert
#include <cctype>               // for isalpha
#include <cstdio>               // for size_t, fgets
#include <cstring>              // for memcpy
#include <exception>            // for exception
#include <iterator>             // for distance
#include <numeric>              // for iota
#include <optional>             // for operator==, optional
#include <set>                  // for operator==
#include <sstream>              // for ostringstream, operator<<, operator==
#include <stdexcept>            // for runtime_error, logic_error
#include <string>               // for string, stoi, stoul, stod, to_string, operator==, operator+
#include <utility>              // for move

using namespace std;
using namespace Decoder;
using namespace Podd;
using OptU_t = std::optional<UInt_t>;
using OptI_t = std::optional<Int_t>;
using OptD_t = std::optional<Double_t>;

namespace {
constexpr UInt_t defaultDT = 4;
}

//_____________________________________________________________________________
ScalerEvtHandler::ScalerEvtHandler( const char* name,
                                          const char* description )
  : THaEvtTypeHandler(name, description)
  , fEvtCount(0)
  , fEvtNum(0)
  , fScalerTree(nullptr)
  , fDeltaT(defaultDT)
{}

//_____________________________________________________________________________
ScalerEvtHandler::~ScalerEvtHandler()
{
  // The tree object is owned by ROOT since it gets associated wth the output
  // file, so DO NOT delete it here.
  if (!TROOT::Initialized())
    delete fScalerTree;
}

//_____________________________________________________________________________
Int_t ScalerEvtHandler::Analyze(THaEvData *evdata)
{
  if( !IsMyEvent(evdata->GetEvType()) && !evdata->IsPhysicsTrigger() )
    return -1;

  Clear();

  // Parse the data, load local data arrays.

  if( Int_t ret = Decode(evdata); ret <= 0 )
    return ret;

  // Copy decoded data to the defined variables
  // Single variables
  for( auto& var: fScalarVars ) {
    assert(var.index < fScalers.size());  // else bug in SetIndices
    var.Fill(fScalers[var.index].get());
  }
  // Arrays
  for( auto& var: fArrayVars ) {
    UInt_t pos = 0;
    for( auto idx: var.idxlist )
      pos = var.Fill(pos, fScalers[idx].get());
  }

  fEvtCount++;
  fEvtNum = evdata->GetEvNum();

  if( fScalerTree )
    fScalerTree->Fill();

  return kOK;
}

//_____________________________________________________________________________
TBranch* ScalerEvtHandler::MakeBranch( const string& name, const THaVar* var ) const
{
  // Build the ROOT tree leaf description: varname/type or varname[size]/type
  string tinfo = name;
  if( var->IsArray()  )
    tinfo += '[' + to_string(var->GetLen() ) + ']';
  auto vartype = var->GetType();
  if( vartype == kUInt || vartype == kUIntP) {
    tinfo += "/i"; // uint32_t
  } else {
    assert(vartype == kFloat || vartype == kFloatP);
    tinfo += "/F"; // float
  }
  // Add a branch for this variable to the tree
  // ROOT really wants a non-const pointer to the data ... pray and hope
  auto* branch = fScalerTree->Branch(
    name.c_str(), const_cast<void*>(var->GetDataPointer()), tinfo.c_str());
  if( !branch )
    Warning("Begin", "Cannot create tree branch \"%s\". "
            "Should never happen. Call expert.", name.c_str());
  return branch;
}

//_____________________________________________________________________________
Int_t ScalerEvtHandler::Begin( THaRunBase* r )
{
  Int_t ret = THaEvtTypeHandler::Begin(r);

  // Add all the scaler variables (scalerloc) to a dedicated scaler tree
  if( !fScalerTree ) {
    TString sname1 = "TS";
    TString sname2 = sname1 + fName;
    TString sname3 = fName + " Scaler Data";

    fScalerTree = new TTree(sname2,sname3);
    fScalerTree->SetAutoSave(200'000'000);

    TString name = "evnum";
    TString tinfo = name + "/l";  // uint64_t
    fScalerTree->Branch(name, &fEvtNum, tinfo, 4000);

    name = "evcount";
    tinfo = name + "/l";          // uint64_t
    fScalerTree->Branch(name, &fEvtCount, tinfo, 4000);

    for( const auto& var: fScalarVars ) {
      assert(var.var); // else bug in DefVars
      assert(&var.count == var.var->GetDataPointer());
      MakeBranch(var.name, var.var);
    }
    for( const auto& arr: fArrayVars ) {
      assert(arr.var); // else bug in DefVars
      assert(arr.pCount == arr.var->GetDataPointer());
      MakeBranch(arr.name, arr.var);
    }
    //TODO chan, slot, crate branches. In separate tree with just one "event"?
  }
  return ret;
}

//_____________________________________________________________________________
void ScalerEvtHandler::Clear( Option_t* string )
{
  THaEvtTypeHandler::Clear(string);

  for( auto& s: fScalers )
    s->Clear();
}

//_____________________________________________________________________________
Int_t ScalerEvtHandler::End( THaRunBase* r )
{
  if( fScalerTree )
    fScalerTree->Write();

  return THaEvtTypeHandler::End(r);
}

//_____________________________________________________________________________
// Helper function for Init()
namespace {
GenScaler* MakeScaler( Int_t model )
{
  // Dynamically create requested decoder module with model number 'model' and
  // ensure that it is a Decoder::GenScaler.
  // Adapted from THaCrateMap::SlotInfo_t::LoadModule.

  // FIXME avoid code duplication -> move this to the decoder

  const char* const here = "THaScalerEvtHandler::MakeScaler";

  auto& moduletypes = Module::fgModuleTypes();
  auto found = moduletypes.find(model);
  if( found == moduletypes.end() ) {
    Error(here, "Decoder module type %d not defined. Crates using this "
          "model will not be decoded. Load the required library.", model);
    return nullptr;
  }
  const auto& modtype = *found;
  assert( modtype.fModel == model );

  // Get the ROOT class for this type, if not already done
  if( !modtype.fTClass ) {
    modtype.fTClass = TClass::GetClass( modtype.fClassName );
    if( !modtype.fTClass ) {
      Error( here, "No ROOT dictionary for class %s. "
             "Coding error. Call expert.", modtype.fClassName );
      return nullptr;
    }
    // Equivalent of dynamic_cast
    const char* const kBaseClassName = "Decoder::GenScaler";
    if( !modtype.fTClass->IsLoaded() ||
        !modtype.fTClass->InheritsFrom(kBaseClassName) ) {
      Error(here, "Class %s does not inherit from %s. "
            "Coding error. Call expert.", modtype.fClassName, kBaseClassName);
      return nullptr;
    }
  }
  assert( modtype.fTClass );

  // If necessary, create new module instance of this type
  auto* module = static_cast<GenScaler*>( modtype.fTClass->New() );
  if( !module ) {
    Error( here, "Failed to make Module %s. Coding error. Call expert.",
           modtype.fClassName );
      return nullptr;
  }

  // Initialize the module
  try {
    module->Init();
  }
  catch( const exception& e ) {
    Error(here, "Failed to initialize module %s: %s",
          modtype.fClassName, e.what());
    delete module;
    return nullptr;
  }
  // Further initialization is done at a later time
  return module;
}
}

//_____________________________________________________________________________
void ScalerEvtHandler::ParseMap( const vector<string>& words )
{
  // Parse cratemap-style map of scaler modules assigned to crate/slot locations
  // to be identified by traditional header word matching (key = "map") or
  // through bank decoding (key = "bank").

  const char* const here = "ParseMap";

  if( words.size() < 5 )
    return;
  const auto& key = words[0];
  UInt_t icrate = stoul(words[1]);
  UInt_t islot  = stoul(words[2]);
  Int_t  imodel = stoi (words[3]);
  OptU_t header, mask;
  OptI_t clkref, bank;

  if( key == "map" ) {
    if( words.size() < 6 ) {
      Error("ParseMap", "\"map\" definition needs at least 6 parameters");
      return;
    }
    header = stoul(words[4], nullptr, 16); // hex
    mask   = stoul(words[5], nullptr, 16); // hex
    if( words.size() >= 7 )
      clkref = stoi(words[6]);
  }
  else if( key == "bank" ) {
    bank   = stoi(words[4]);
    if( words.size() >= 6 )
      clkref = stoi(words[5]);
  } else {
    assert(false);  // else error in ReadDatabase
#ifdef NDEBUG
    Error( "ParseMap", "Unknown keyword %s", key.c_str() );
    return;
#endif
  }

  // Create the requested decoder module
  unique_ptr<GenScaler> scaler{MakeScaler(imodel)};
  // Sanity checks
  if( !scaler ) {
    ostringstream oss;
    oss << "Unable to create GenScaler module type " << imodel
      << ". Call expert.";
    throw runtime_error(oss.str()); // probably didn't load relevant library
  }
  if( scaler->GetModelNum() != imodel ) {
    ostringstream oss;
    oss << "GenScaler module type " << imodel
      << " model number mismatch. Coding error. Call expert.";
    throw logic_error(oss.str());  // probably bug in Module::Init()
  }
  // Check for duplicates - each crate/slot can only be assigned one module
  if( FindScaler(icrate, islot) != fScalers.end() ) {
    Error(here, "Duplicate scaler definition \"%s %u %u %d\". "
          "Crate/slot must be unique. Ignoring line.",
          words[0].c_str(), icrate, islot, imodel);
    return;
  }
  // Configure the module
  scaler->SetSlot(icrate, islot );
  scaler->SetModel(imodel);
  if( header && mask ) {
    scaler->SetHeader(header.value(), mask.value());
  }
  else if( bank ) {
    // SetBank() sets header and mask, too
    scaler->SetBank(bank.value());
  }
  // Hijack unused mode variable for the clock reference
  scaler->SetMode(clkref.value_or(-1));
  fScalers.push_back(std::move(scaler));
}

//_____________________________________________________________________________
void ScalerEvtHandler::ParseClock( const vector<string>& words )
{
  if( words.size() < 2 )
    return; // Quietly ignore bare "clock" lines
  OptU_t icrate, islot, ichan;
  OptI_t clkref;
  OptD_t deltaT, clkfreq;
  if( words.size() >= 6 ) {
    clkref  = stoi(words[1]);
    icrate  = stoul(words[2]);
    islot   = stoul(words[3]);
    ichan   = stoul(words[4]);
    clkfreq = stod(words[5]);
  } else if ( words.size() == 2 ) {
    deltaT  = stod(words[1]);
  } else {
    throw runtime_error("Garbled clock definition. Check database.");
  }
  if( deltaT )
    fDeltaT = deltaT.value();
  else {
    auto iclk = clkref.value();
    auto dup = ranges::find_if(fClocks, [iclk](const auto& c)
      { return iclk == c.iref; });
    if( dup != fClocks.end() ) {
      ostringstream oss;
      oss << "Duplicate clock number \"" << words[0] << " " << iclk << " "
        << icrate.value() << " " << islot.value() << " " << ichan.value()
        << clkfreq.value() << "\". Ignoring line.";
      throw runtime_error(oss.str());
    }
    fClocks.push_back({
      .iref = clkref.value(), .icrate = icrate.value(), .islot = islot.value(),
      .ichan = ichan.value(), .freq = clkfreq.value()
    });
  }
}

//_____________________________________________________________________________
namespace {
Int_t ParseVariableLine( const vector<string>& words, UInt_t nameidx,
  OptU_t& icrate, OptU_t& islot, OptU_t& ichan, OptU_t& imode, OptI_t& ibank )
{
  bool is_variable = false;
  const auto& key = words[0];
  if( key == "variable" ) {
    if( nameidx < 4 || nameidx > 5 ) {
      Warning("ParseVariable", "\"variable\" definition needs "
              "3 or 4 numbers, followed by the name. Variable not defined.");
      return -1;
    }
    icrate = stoul(words[1]);
    islot  = stoul(words[2]);
    ichan  = stoul(words[3]);
    if( nameidx == 5 )
      imode  = stoul(words[4]);
    is_variable = true;
  } else if( key == "array" ) {
    if( nameidx > 4 ) {
      Warning("ParseVariable", "Malformed \"array\" definition. "
              "There must be at most 3 numbers, followed by the name. "
              "Variable not defined.");
      return -1;
    }
    switch( nameidx ) {
      case 4: // 3 numbers: crate, slot, mode
        islot = stoul(words[2]);
        [[fallthrough]];
      case 3: // 2 numbers crate, mode
        icrate = stoul(words[1]);
        [[fallthrough]];
      case 2: // 1 number: mode
        imode = stoul(words[nameidx - 1]);
        [[fallthrough]];
      default: // no numbers: ok
        break;
    }
  } else if( key == "bankarray" ) {
    if( nameidx == 1 || nameidx > 5 ) {
      Warning("ParseVariable", "Malformed \"bankarray\" definition. "
              "There must between 1 and 4 numbers, followed by the name. "
              "Variable not defined.");
      return -1;
    }
    switch( nameidx ) {
      case 5: // 4 numbers: bank, crate, slot, mode
        islot = stoul(words[3]);
        [[fallthrough]];
      case 4: // 3 numbers bank, crate, mode
        icrate = stoul(words[2]);
        [[fallthrough]];
      case 3: // 2 number: bank, mode
        imode = stoul(words[nameidx - 1]);
        [[fallthrough]];
      case 2: // 1 number: bank
        ibank = stoul(words[1]);
        break;
      default: // no numbers: disallowed, should have been caught above
        assert(false);
#ifdef NDEBUG
        return -1;
#endif
    }
  } else {
    assert(false); // else error in ReadDatabase
#ifdef NDEBUG
    Error("ParseVariable", "Unknown keyword %s", key.c_str());
    return -1;
#endif
  }
  return is_variable;
}

//_____________________________________________________________________________
UInt_t FindNameIdx( const vector<string>& words )
{
  UInt_t nameidx;
  for( nameidx = 1; nameidx < 6 && nameidx < words.size() &&
                  !isalpha(words[nameidx][0]); ++nameidx ) {}
  if( nameidx == words.size() ) {
    Error("ParseVariable", "No variable name given. Check database.");
    return 0;
  }
  return nameidx;
}

//_____________________________________________________________________________
ScalerEvtHandler::EPick GetPick( OptU_t icrate, OptU_t islot )
{
  using enum ScalerEvtHandler::EPick;
  ScalerEvtHandler::EPick ipick;
  if( islot ) {
    ipick = kSlot;
    assert(icrate);
  } else if( icrate ) {
    ipick = kCrate;
  } else
    ipick = kAll;
  return ipick;
}
}

//_____________________________________________________________________________
void ScalerEvtHandler::ParseVariable( const vector<string>& words )
{
  if( words.size() < 2 )
    return;

  // Find the variable name, which must start with a letter
  UInt_t nameidx = FindNameIdx(words);
  if( nameidx == 0 )
    return;

  OptU_t icrate, islot, ichan, imode;
  OptI_t ibank;
  Int_t ret = ParseVariableLine(words, nameidx, icrate, islot, ichan,
                                imode, ibank);
  if( ret < 0 )
    return;
  bool is_variable = ret > 0;

  string varname = words[nameidx];
  // The leftover words on the line are taken as the variable description,
  // prefixed by the name of this event handler.
  string vardesc(GetName());
  for( size_t j = nameidx + 1; j < words.size(); j++ )
    vardesc += " " + words[j];
  // Determine whether this variable reports raw counts, the rate, or both.
  // ".count" and ".rate" will be appended to the variable name, as applicable.
  bool count_mode = imode && (imode.value() == 0 || imode.value() == 1);
  // if not set, assume "rate"
  bool rate_mode  = !imode || imode.value() == 0 || imode.value() == 2;
  if( !count_mode && !rate_mode ) {
    Error("ParseVariable", "Unrecognized data type %u. Must be 1 (counts),"
          "2 (rate), or 0 (both). Fix database.", imode.value());
    return;
  }

  if( is_variable ) {
    assert(icrate && islot && ichan);
    if( count_mode )
      fScalarVars.emplace_back(varname, vardesc, icrate.value(), islot.value(),
        ichan.value(), kCount);
    if( rate_mode )
      fScalarVars.emplace_back(varname, vardesc, icrate.value(), islot.value(),
        ichan.value(), kRate);

  } else {
    auto ipick = GetPick(icrate, islot);
    if( count_mode )
      fArrayVars.emplace_back(varname, vardesc, icrate.value_or(0),
        islot.value_or(0), kCount, ibank.value_or(-1), ipick);
    if( rate_mode )
      fArrayVars.emplace_back(varname, vardesc, icrate.value_or(0),
        islot.value_or(0), kRate, ibank.value_or(-1), ipick);
  }
}

//_____________________________________________________________________________
Int_t ScalerEvtHandler::AssignNormScaler()
{
  const char* const here = "AssignNormScaler";

  //TODO allow running without a clock, just estimate rates with deltaT
  if( fClocks.empty() ) {
    Error(here, "No clocks defined. Fix database.");
    return kInitError;
  }
  bool multi_clocks = fClocks.size() > 1;
  // Connect scalers to clock(s)
  for( auto& scaler : fScalers ) {
    UInt_t clkidx;
    Int_t clkref = scaler->GetMode();
    if( clkref < 0 ) {
      if( multi_clocks ) {
        Error(here, "Scaler in crate/slot %u/%u defined to use "
              "default clock, but multiple clocks defined. Fix database.",
              scaler->GetCrate(), scaler->GetSlot());
        return kInitError;
      }
      clkidx = 0;
    } else {
      auto iclock = ranges::find_if(fClocks, [clkref]( const auto& c )
        { return c.iref == clkref; });
      if( iclock == fClocks.end() ) {
        Error(here, "Clock reference %d for scaler in crate/slot %u/%u does "
              "not exist. Fix database.",
              clkref, scaler->GetCrate(), scaler->GetSlot());
        return kInitError;
      }
      clkidx = std::distance(fClocks.begin(), iclock);
    }
    assert(clkidx < fClocks.size());
    const auto& [iref, icrate, islot, ichan, freq] = fClocks[clkidx];
    assert(iref == clkref);
    const auto iclkscal = FindScaler(icrate, islot);
    if( iclkscal == fScalers.end() ) {
      Error(here, "No scaler for clock definition \"%s %d %u %u %u %lf\". "
              "Fix database.", "clock", iref, icrate, islot, ichan, freq);
      return kInitError;
    }
    // At this point, the current scaler refers to an existing clock definition,
    // and that clock definition references an existing scaler
    if( icrate == scaler->GetCrate() && islot == scaler->GetSlot() ) {
      // Self-reference
      scaler->SetClock(fDeltaT, ichan, freq);
    } else {
      scaler->LoadNormScaler(iclkscal->get());
    }
  }
  return kOK;
}

//_____________________________________________________________________________
Int_t ScalerEvtHandler::Decode( THaEvData* evdata )
{
  // For each scaler, get its crate/slot. Query the crate map whether it is
  // a bank-type or a header-type slot. Decode accordingly

  const char* const here = "Decode";

  bool isCoda3 = evdata->GetDataVersion() > 2;
  auto* codaevent = dynamic_cast<CodaDecoder*>(evdata);
  if( !codaevent && isCoda3 ) {
    // For CODA3, we need the CodaDecoder to handle bank data
    Error(here, "Decoder is not a CodaDecoder. Cannot proceed.");
    return -1;
  }
  const auto* crmap = evdata->GetCrateMap();
  if( !crmap ) {
    Error(here, "Crate map not available. Cannot proceed.");
    return -2;
  }

  bool found = false;
  // TODO Use lists of used banks, used crates
  for( auto& scaler: fScalers ) {
    bool do_bank = false;
    if( isCoda3 ) {
      if( crmap->isAllBanks(scaler->GetCrate()) )
        do_bank = true;
    }
    Int_t ret;
    if( do_bank ) {
      ret = DecodeBank(codaevent, scaler.get());
    } else {
      ret = DecodeRoc(evdata, scaler.get());
    }
    if( ret > 0 )
      found = true;
  }
  return found;
}

//_____________________________________________________________________________
Int_t ScalerEvtHandler::DecodeBank( CodaDecoder* codaevent,
                                       GenScaler* scaler )
{
  // Decode given scaler whose data are bank-structured
  // Return 1 if bank found and scaler has data, 0 otherwise

  const char* whereami = Here("DecodeBank"); //TODO
  using enum CodaDecoder::BankInfo::EDataType;

  // Unpack CODA 3 bank structure
  auto crate = scaler->GetCrate();
  auto bank = scaler->GetBank();
  assert(bank > 0);
  auto bankloc = codaevent->FindBank(crate, bank);
  if( bankloc.pos == 0 )
    // Bank not found. This is expected. Not all events have all banks.
    return 0;

  // Find the payload for this crate/bank
  auto bankinfo = codaevent->GetBank(bankloc);

  // Check for errors
  if( !bankinfo ) {
    ostringstream ostr;
    ostr << whereami << ": CODA3 bank decoding error \""
      << bankinfo.Errtxt() << "\"";
    throw CodaDecoder::coda_format_error(ostr.str());
  }
  auto type = bankinfo.GetDataType();
  if( type != kUInt && type != kInt ) {
    ostringstream ostr;
    ostr << whereami << ": CODA3 bank does not contain 32-bit integer "
      "data, found " << bankinfo.Typtxt();
    throw CodaDecoder::coda_format_error(ostr.str());
  }

  // Let the scaler module pull its slot data out of the bank
  scaler->LoadSlot(nullptr, codaevent->GetRawDataBuffer(),
                   bankinfo.pos_, bankinfo.len_);

  return scaler->IsDecoded();
}

//_____________________________________________________________________________
Int_t ScalerEvtHandler::DecodeRoc( THaEvData* evdata, GenScaler* scaler )
{
  // Decode given scaler from legacy ROC-structured data
  // Return 1 if data found and scaler has data, 0 otherwise

  const UInt_t *p, *pstop;
  if( evdata->IsPhysicsTrigger() ) {
    // For physics events, scan the ROC block
    auto crate = scaler->GetCrate();
    UInt_t roclen = evdata->GetRocLength(crate);
    if( roclen == 0 )
      return 0;
    p = evdata->GetRawDataBuffer(crate) + 2;
    pstop = p + roclen;
  } else {
    // User events (presumably) do not have a ROC structure. Scan the
    // entire event buffer (as it is done in Bob's old code).
    // FIXME this does not work with scalers in multiple crates.
    // FIXME so ensure that only one crate specified in map
    p = evdata->GetRawDataBuffer();
    pstop = p + evdata->GetEvLength() - 1; // last word of data
  }
  scaler->LoadSlot(nullptr, p, pstop);

  return scaler->IsDecoded();
}

//_____________________________________________________________________________
void ScalerEvtHandler::SetIndices()
{
  // Associate ScalerVariable objects with single scaler modules
  for( auto& svar: fScalarVars ) {
    auto found = ranges::find_if(fScalers,
      [&]( const auto& s ) { return svar.Match(s.get()); });
    bool good = found != fScalers.end();
    if( good ) {
      svar.index = std::distance(fScalers.begin(), found);
      auto nchan = found->get()->GetNumChan();
      if( svar.ichan >= nchan ) {
        Warning("ParseVariable", "Requested channel %u for scaler in crate = "
                "%u, slot = %u out of range (max %u). Variable \"%s\" not "
                "defined.", svar.ichan, svar.icrate, svar.islot, nchan,
                svar.name.c_str());
        good = false;
      }
    } else {
      Warning("ParseVariable", "No scaler module in crate = %u, "
              "slot = %u available. Variable \"%s\" not defined.",
              svar.icrate, svar.islot, svar.name.c_str());
    }
    if( !good )
      svar.name.clear();  // simple way to flag invalid definitions
  }
  erase_if(fScalarVars, [](const auto& svar) { return svar.name.empty(); });

  // ScalerArray variables
  for( auto& sarr: fArrayVars ) {
    // Count total number of channels
    UInt_t totchan = 0;
    bool matched = false;
    auto& idxs = sarr.idxlist;
    idxs.clear();
    idxs.reserve(fScalers.size());
    for( UInt_t i = 0; i < fScalers.size(); ++i ) {
      const GenScaler* s = fScalers[i].get();
      if( sarr.Match(s) ) {
        totchan += s->GetNumChan();
        idxs.push_back(i);
        matched = true;
      }
    }
    if( !matched || totchan == 0 ) {
      Warning("ParseVariable", "No scaler module matches array variable "
              "\"%s\". Variable not defined.", sarr.name.c_str());
      sarr.name.clear();
      continue;
    }
    // Allocate memory for data and indices, as needed
    sarr.Allocate(totchan);
    // Fill index array(s)
    UInt_t pos = 0;
    for( auto idx: idxs )
      pos = sarr.FillIndices(pos, fScalers[idx].get());
  }
  erase_if(fArrayVars, [](const auto& sarr) { return sarr.name.empty(); });
}

//_____________________________________________________________________________
Int_t ScalerEvtHandler::ReadDatabase( const TDatime& date )
{
  // Parse the map file which defines what scalers exist and the global variables.
  CFile file = OpenFile(date);
  if( !file )
    return kFileError;

  constexpr int LEN = 256;
  char cbuf[LEN];

  while( fgets(cbuf, LEN, file) ) {
    string line{cbuf};
    auto pos = line.find('#');
    if( pos != string::npos )
      line.erase(pos);
    Podd::Trim(line);  // also erases trailing newline that fgets keeps
    if( line.empty() )
      continue;

    const vector<string> words = vsplit(line);
    try {
      if( words.front() == "map" || words.front() == "bank" ) {
        ParseMap(words);
      } else if( words.front() == "clock" ) {
        ParseClock(words);
      } else if( words.front() == "variable" || words.front() == "array" ||
                 words.front() == "bankarray" ) {
        ParseVariable(words);
      }
    } catch( const runtime_error& e ) {
      Error(Here("ReadDatabase"), "%s\n\"%s\"", e.what(), line.c_str());
      return kInitError;
    } catch( const std::exception& e ) {
      Error(Here("ReadDatabase"), "Caught %s parsing database line\n"
            "\"%s\"", e.what(), line.c_str());
      return kInitError;
    }
  }
  return kOK;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus ScalerEvtHandler::Init(const TDatime& date)
{
  fStatus = THaEvtTypeHandler::Init(date);  // calls ReadDatabase
  if( fStatus != kOK )
    return fStatus;

  // Set this from the analysis script
  //  AddEvtType(140);  // what events to look for

  // Identify indices of scalers[] vector to variables.
  SetIndices();

  // Call LoadNormScaler or SetClock after scalers created
  Int_t ret = AssignNormScaler();
  if( ret != kOK )
    return fStatus = kInitError;

  // Define global variables on all scalerloc objects
  fStatus = DefVars();

  return fStatus;
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus ScalerEvtHandler::DefVars()
{
  // Called after ParseVariables and SetIndices

  if( !gHaVars ) {
    Error("DefVars", "No gHaVars ?!  Well, that's a problem !!");
    return kInitError;
  }
  size_t Nvars = fScalarVars.size() + fArrayVars.size();
  if( Nvars == 0 )
    return kOK;

  // Variables
  for( auto& svar: fScalarVars ) {
    // As usual, prefix global variables with the module's prefix
    string varname = GetPrefix() + svar.name;
    switch(svar.ikind) {
      case kCount:
        svar.var = gHaVars->Define(varname.c_str(), svar.description.c_str(),
                                   svar.count);
        break;
      case kRate:
        svar.var = gHaVars->Define(varname.c_str(), svar.description.c_str(),
                                   svar.rate);
        break;
    }
  }

  // Arrays
  for( auto& sarr: fArrayVars) {
    assert(sarr.size > 0 );   // else bug in SetIndices()
    // The size of the array variables is fixed; it is set once at Init() time
    string varname = GetPrefix() + sarr.name;
    string subscript = '[' + to_string(sarr.size) + ']';
    string arrname = varname + subscript;
    switch( sarr.ikind ) {
      case kCount:
        sarr.var = gHaVars->Define(arrname.c_str(), sarr.description.c_str(),
                                   sarr.pCount);
        break;
      case kRate:
        sarr.var = gHaVars->Define(arrname.c_str(), sarr.description.c_str(),
                                   sarr.pRate);
        break;
    }
    // Define chan, slot, crate variables. These currently hold constant data
    // that never change from event to event, which wastes a lot of space, but
    // they are the most convenient way to recover the source of the scaler
    // info from the parallel .count and/or .rate arrays.
    // On the other hand, ROOT file compression will be very effective on
    // corresponding branches.
    // TODO can these arrays be written just once in some kind of header?

    // Drop trailing ".count" or ".rate" in name
    ChopPrefix( varname);
    // Similarly, drop trailing "(count)" or "(rate)" in description
    string vardesc{sarr.description};
    auto pos = vardesc.rfind('(');
    if( pos != string::npos )
      vardesc.erase(pos);
    string auxname, auxdesc;
    switch( sarr.ipick ) {
      case kAll:
        assert(sarr.pCrate);
        auxname = varname + "crate"; auxname += subscript;
        auxdesc = vardesc + "(crate number)";
        gHaVars->Define(auxname.c_str(), auxdesc.c_str(), sarr.pCrate);
        [[fallthrough]];
      case kCrate:
        assert(sarr.pSlot);
        auxname = varname + "slot"; auxname += subscript;
        auxdesc = vardesc + "(slot number)";
        gHaVars->Define(auxname.c_str(), auxdesc.c_str(), sarr.pSlot);
        [[fallthrough]];
      case kSlot:
        assert(sarr.pChan);
        auxname = varname + "chan"; auxname += subscript;
        auxdesc = vardesc + "(channel number)";
        gHaVars->Define(auxname.c_str(), auxdesc.c_str(), sarr.pChan);
        break;
    }
  }
  return kOK;
}

//_____________________________________________________________________________
decltype(ScalerEvtHandler::fScalers)::iterator
ScalerEvtHandler::FindScaler( UInt_t icrate, UInt_t islot )
{
  return ranges::find_if(fScalers, [icrate,islot]( const auto& s ) {
    return s->GetCrate() == icrate && s->GetSlot() == islot;
  });
}

//_____________________________________________________________________________
string ScalerEvtHandler::Variable::KindName( EKind kind )
{
  constexpr const char* kCountName = ".count";
  constexpr const char* kRateName = ".rate";
  return kind == kCount ? kCountName : kRateName;
}

//_____________________________________________________________________________
std::string ScalerEvtHandler::Variable::KindDesc( EKind kind )
{
  constexpr const char* kCountName = " (count)";
  constexpr const char* kRateName = " (rate)";
  return kind == kCount ? kCountName : kRateName;
}

//_____________________________________________________________________________
void ScalerEvtHandler::ScalarVariable::Fill( const GenScaler* scaler )
{
  assert(ichan < scaler->GetNumChan()); // else bug in SetIndices
  switch( ikind ) {
    case kCount:
      count = scaler->GetData(ichan);
      break;
    case kRate:
      rate = static_cast<Float_t>(scaler->GetRate(ichan));
      break;
  }
}

//_____________________________________________________________________________
ScalerEvtHandler::ArrayVariable::ArrayVariable(
      const std::string& nm, const std::string& desc, UInt_t cr, UInt_t sl,
      EKind kind, Int_t bank, EPick pick )
  : Variable( nm, desc, cr, sl, kind )
  , size(0)  // set later
  , ibank(bank)
  , ipick(pick)
  , pCount(nullptr)
  , pChan(nullptr)
  , pSlot(nullptr)
  , pCrate(nullptr)
{}

// We have to implement these move operations to be able to use ArrayVariable
// objects in a std::vector
//_____________________________________________________________________________
ScalerEvtHandler::ArrayVariable::ArrayVariable( ArrayVariable&& rhs ) noexcept
  : Variable(std::move(static_cast<Variable&&>(rhs))) // intentionally slices
  , size(rhs.size)
  , ibank(rhs.ibank)
  , ipick(rhs.ipick)
  , idxlist(std::move(rhs.idxlist))
  , pCount(rhs.pCount)
  , pChan(rhs.pChan)
  , pSlot(rhs.pSlot)
  , pCrate(rhs.pCrate)
{
  rhs.pCount = nullptr;
  rhs.pChan = nullptr;
  rhs.pSlot = nullptr;
  rhs.pCrate = nullptr;
}

//_____________________________________________________________________________
ScalerEvtHandler::ArrayVariable&
  ScalerEvtHandler::ArrayVariable::operator=( ArrayVariable&& rhs ) noexcept
{
  if( this == &rhs )
    return *this;
  Variable::operator=(static_cast<Variable&&>(rhs)); // intentionally slices
  size = rhs.size;
  ibank = rhs.ibank;
  ipick = rhs.ipick;
  idxlist = std::move(rhs.idxlist);
  pCount = rhs.pCount;
  pChan = rhs.pChan;
  pSlot = rhs.pSlot;
  pCrate = rhs.pCrate;
  rhs.pCount = nullptr;
  rhs.pChan = nullptr;
  rhs.pSlot = nullptr;
  rhs.pCrate = nullptr;
  return *this;
}

//_____________________________________________________________________________
bool ScalerEvtHandler::ArrayVariable::Match( const GenScaler* scaler ) const
{
  // Check if this array variable is associated with the given scaler

  if( MatchesAll() )
    return true;

  switch( ipick ) {
    case kSlot:
      if( islot != scaler->GetSlot() )
        return false;
      [[fallthrough]];
    case kCrate:
      if( icrate != scaler->GetCrate() )
        return false;
      [[fallthrough]];
    case kAll:
      if( ibank > 0 && ibank != scaler->GetBank() )
        return false;
      break;
  }
  return true;
}

//_____________________________________________________________________________
void ScalerEvtHandler::ArrayVariable::Allocate( const UInt_t numchan )
{
  // Allocate raw C-arrays for use with ROOT tree branches

  Deallocate(); // avoid memory leaks
  size = numchan;
  if( ikind == kCount ) {
    pCount = new UInt_t[numchan];
  } else {
    assert(ikind == kRate);
    pRate = new Float_t[numchan];
  }
  // Don't allocate redundant index arrays that would hold all identical values
  switch( ipick ) {
    case kAll:
      pCrate = new UShort_t[numchan];
      [[fallthrough]];
    case kCrate:
      pSlot = new Byte_t[numchan];
      [[fallthrough]];
    case kSlot:
      pChan = new Byte_t[numchan];
      break;
  }
}

//_____________________________________________________________________________
void ScalerEvtHandler::ArrayVariable::Deallocate()
{
  // Deallocate the C-arrays that ROOT so loves

  // This is not strictly necessary, I think - UInt_t and Float_t have the
  // same size, so either deletion should be equivalent. But who knows ...
  if( ikind == kCount )
    delete [] pCount;
  else {
    assert(ikind == kRate);
    delete [] pRate;
  }
  pCount = nullptr;
  delete [] pCrate; pCrate = nullptr;
  delete [] pSlot;  pSlot = nullptr;
  delete [] pChan;  pChan = nullptr;
}

//_____________________________________________________________________________
UInt_t ScalerEvtHandler::ArrayVariable::Fill( UInt_t pos,
  const GenScaler* scaler ) const
{
  UInt_t nchan = scaler->GetNumChan();
  assert(pos + nchan <= size);
  assert(pChan[pos] == 0);
  assert(!pSlot || pSlot[pos] == scaler->GetSlot());
  assert(!pCrate || pCrate[pos] == scaler->GetCrate());
  if( ikind == kCount ) {
    const auto& counts = scaler->GetCounts();
    if( counts.size() != nchan )
      throw logic_error("Unexpected scaler count array size mismatch");
    memcpy(pCount + pos, counts.data(), nchan * sizeof(UInt_t));
  } else {
    assert(ikind == kRate);
    const auto& rates = scaler->GetRates();
    if( rates.size() != nchan )
      throw logic_error("Unexpected scaler rates array size mismatch");
    for( UInt_t i = 0; i < nchan; i++ )
      pRate[pos + i] = static_cast<Float_t>(rates[i]);
  }
  return pos + nchan;
}

//_____________________________________________________________________________
UInt_t ScalerEvtHandler::ArrayVariable::FillIndices( UInt_t pos,
  const GenScaler* scaler ) const
{
  UInt_t nchan = scaler->GetNumChan();
  assert(pos + nchan <= size);
  switch( ipick ) {
    case kAll:
      std::fill_n(pCrate + pos, nchan, scaler->GetCrate());
      [[fallthrough]];
    case kCrate:
      std::fill_n(pSlot + pos, nchan, scaler->GetSlot());
      [[fallthrough]];
    case kSlot:
      std::iota(pChan + pos, pChan + pos + nchan, 0);
      break;
  }
  return pos + nchan;
}

//_____________________________________________________________________________
#if ROOT_VERSION_CODE < ROOT_VERSION(6,36,0)
ClassImp(THaScalerEvtHandler)
#endif
