/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// -------------------------------------------------------------------------
// -----                      R3BTTTXMappedData source file            -----
// -------------------------------------------------------------------------

#include "R3BTTTXMappedData.h"
#include <fmt/core.h>

R3BTTTXMappedData::R3BTTTXMappedData(uint8_t DetID, uint8_t StripID, int32_t time, int32_t energy, bool pu, bool ov)
    : fDetID(DetID)
    , fStripID(StripID)
    , fTime(time)
    , fEnergy(energy)
    , fPileup(pu)
    , fOverflow(ov)
{
}

std::string R3BTTTXMappedData::toString() const
{
    return fmt::format("DetID: {}, StripID: {}, Time: {}, Energy: {}, Pileup: {}, Overflow: {}",
                       GetDetID(),
                       GetStripID(),
                       GetTime(),
                       GetEnergy(),
                       GetPileupStatus(),
                       GetOverflowStatus());
}

void R3BTTTXMappedData::Print(const Option_t*) const { std::cout << *this << std::endl; }

std::ostream& operator<<(std::ostream& os, const R3BTTTXMappedData& data)
{
    os << data.toString();
    return os;
}

ClassImp(R3BTTTXMappedData)
