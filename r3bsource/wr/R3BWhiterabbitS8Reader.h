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

#ifndef R3BWhiterabbitS8Reader_H
#define R3BWhiterabbitS8Reader_H

#include "R3BReader.h"
#include <Rtypes.h>

struct EXT_STR_h101_WRS8_t;
typedef struct EXT_STR_h101_WRS8_t EXT_STR_h101_WRS8;

class TClonesArray;
class R3BEventHeader;

class R3BWhiterabbitS8Reader : public R3BReader
{
  public:
    // Standard constructor
    R3BWhiterabbitS8Reader(EXT_STR_h101_WRS8*, size_t, UInt_t);

    // Destructor
    virtual ~R3BWhiterabbitS8Reader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t Read() override;

    // Reset
    virtual void Reset() override;

    /** Accessor to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    // An event counter
    UInt_t fNEvent;
    // Reader specific data structure from ucesb
    EXT_STR_h101_WRS8* fData;
    // Offset of detector specific data in full data structure
    size_t fOffset;
    // The whiterabbit subsystem ID
    UInt_t fWhiterabbitId;
    // A pointer to the R3BEventHeader structure
    R3BEventHeader* fEventHeader;
    // Don't store data for online
    Bool_t fOnline;
    // Output array
    TClonesArray* fArray;

  public:
    ClassDefOverride(R3BWhiterabbitS8Reader, 0);
};
#endif // R3BWhiterabbitS8Reader_H
