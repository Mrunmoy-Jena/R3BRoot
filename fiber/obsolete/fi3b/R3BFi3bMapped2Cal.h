/******************************************************************************
 *   Copyright (C) 2023 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2023 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BFI3BMAPPED2CAL
#define R3BFI3BMAPPED2CAL

#include "R3BBunchedFiberMapped2Cal.h"

class R3BFi3bMapped2Cal : public R3BBunchedFiberMapped2Cal
{
  public:
    R3BFi3bMapped2Cal(Int_t, enum R3BTCalEngine::CTDCVariant);
    R3BFi3bMapped2Cal(Int_t, enum Electronics, enum R3BTCalEngine::CTDCVariant);
    virtual ~R3BFi3bMapped2Cal();

    ClassDef(R3BFi3bMapped2Cal, 1)
};

#endif
