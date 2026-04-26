#ifndef Podd_THaVarList_h_
#define Podd_THaVarList_h_

//////////////////////////////////////////////////////////////////////////
//
// THaVarList
//
//////////////////////////////////////////////////////////////////////////

#include "THashList.h"  // for THashList
#include "THaVar.h"     // for THaVar (ptr only)  // IWYU pragma: keep
#include "VarDef.h"     // for VarDef, RVarDef, etc.
#include <vector>       // for vector
class TString;

class THaVarList : public THashList {
  
public:
  THaVarList();

  // Define() with reference to variable
  template<typename T>
  THaVar* Define( const char* name, const char* descript, const T& var,
                  const Int_t* count = nullptr )
  {
    auto type = Vars::FindType(typeid(T));
    return DefineByType(name, descript, &var, type, count);
  }

  template<typename T>
  THaVar* Define( const char* name, const T& var, const Int_t* count = nullptr )
  { return Define(name, name, var, count); }

  // Define() with pointer to variable
  template<typename T>
  THaVar* Define( const char* name, const char* descript, const T* const& var,
                  const Int_t* count = nullptr )
  {
    auto type = Vars::FindType(typeid(T));
    return DefineByType(name, descript, &var, type, count);
  }

  template<typename T>
  THaVar* Define( const char* name, const T* const& var,
                  const Int_t* count = nullptr )
  { return Define(name, name, var, count); }

  //Avoid ambiguity - cannot specify variable length char array like this
  //Must use the form with description
  THaVar* Define( const char* name, const Char_t* const& var )
  {
    return Define(name, name, var);
  }

  virtual THaVar*  DefineByType( const char* name, const char* desc,
				 const void* loc, VarType type,
				 const Int_t* count,
				 const char* errloc = "DefineByType" );
  virtual THaVar*  DefineByRTTI( const TString& name, const TString& desc,
				 const TString& def, const void* obj,
				 TClass* cl,
				 const char* errloc = "DefineByRTTI" );

  struct DefineVariablesOpts {
    const char* prefix        {""};
    const char* caller        {""};
    const char* def_prefix    {""};
    const char* comment_subst {""};
  };
  virtual Int_t    DefineVariables( const std::vector<VarDef>& list,
                                    DefineVariablesOpts opts );
  virtual Int_t    DefineVariables( const std::vector<RVarDef>& list,
                                    const TObject* obj,
                                    DefineVariablesOpts opts );
  // Legacy methods for backward compatibility
  [[deprecated("Use overload with vector<VarDef>")]]
  virtual Int_t    DefineVariables( const VarDef* list,
				    const char* prefix="",
				    const char* caller="",
                                    const char* def_prefix = "",
                                    const char* comment_subst = "");
  [[deprecated("Use overload with vector<RVarDef>")]]
  virtual Int_t    DefineVariables( const RVarDef* list,
                                    const TObject* obj,
                                    const char* prefix = "",
                                    const char* caller = "",
                                    const char* def_prefix = "",
                                    const char* comment_subst = "");
  virtual THaVar*  Find( const char* name ) const;
  virtual void     PrintFull(Option_t *opt="") const;
  virtual Int_t    RemoveName( const char* name );
  virtual Int_t    RemoveRegexp( const char* expr, Bool_t wildcard = true );

protected:

  ClassDef(THaVarList,2)   //List of analyzer global variables
};

#endif

