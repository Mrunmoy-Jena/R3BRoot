/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2023 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BSyncCheckReader_H
#define R3BSyncCheckReader_H 1

#include "R3BReader.h"
#include <Rtypes.h>

struct EXT_STR_h101_SYNC_CHECK_t;
typedef struct EXT_STR_h101_SYNC_CHECK_t EXT_STR_h101_SYNC_CHECK;
class ext_data_struct_info;
class TClonesArray;

class R3BSyncCheckReader : public R3BReader
{
  public:
    // Standard constructor
    R3BSyncCheckReader(EXT_STR_h101_SYNC_CHECK*, size_t);

    // Destructor
    virtual ~R3BSyncCheckReader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t R3BRead() override;

    // Reset
    virtual void Reset() override;

  private:
    // An event counter
    UInt_t fNEvent;
    // Reader specific data structure from ucesb
    EXT_STR_h101_SYNC_CHECK* fData;
    // Offset of detector specific data in full data structure
    size_t fOffset;
    // Output array
    TClonesArray* fArray;

  public:
    ClassDefOverride(R3BSyncCheckReader, 0);
};
#endif // R3BSyncCheckReader_H