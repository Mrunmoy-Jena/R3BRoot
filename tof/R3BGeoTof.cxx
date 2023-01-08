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

/////////////////////////////////////////////////////////////
// R3BGeoTof
//
// Class for geometry of R3BCAL
//
/////////////////////////////////////////////////////////////
#include "R3BGeoTof.h"
#include "FairGeoNode.h"
#include "FairLogger.h"
#include <iostream>

ClassImp(R3BGeoTof)

    R3BGeoTof::R3BGeoTof()
{
    // Constructor
    fName = "sts";
    maxSectors = 0;
    maxModules = 99;
}

const char* R3BGeoTof::getModuleName(Int_t m)
{
    // Returns the module name of sts number m
    if (m < 0)
    {
        LOG(error) << "R3BGeoTof::getModuleName:: Module number " << m << " not known!";
        return "";
    }
    if (m < 9)
        sprintf(modName, "calstation0%i", m + 1);
    else
        sprintf(modName, "calstation%i", m + 1);
    return modName;
}

const char* R3BGeoTof::getEleName(Int_t m)
{
    // Returns the element name of sts number m
    sprintf(eleName, "cal%i", m + 1);
    return eleName;
}
