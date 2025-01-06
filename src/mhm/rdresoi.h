//
// rdresoi.h
//
//    Copyright (C) 2023 by Wuhan University
//
//    This program belongs to PRIDE PPP-AR which is an open source software:
//    you can redistribute it and/or modify it under the terms of the GNU
//    General Public License (version 3) as published by the Free Software
//Foundation.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License (version 3) for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contributor: Honghai Zhang
//
//
//
// purpose   : read the main body information of the residual file
//
#ifndef RDRESOH_H
#define RDRESOH_H

#include "rdresoh.h"
#define MAX_FILENAME_LENGTH 128
#define MAX_LINE_LENGTH 1024
#define MAX_SATELLITES 10000
#define MAX_SATELLITE_NAME_LENGTH 4

struct SatelliteData {
    float timeInSeconds;
    double phaseResidue;
    double phaseResidue_filtered;
    double pseudoRangeResidue;
    double pseudoRangeResidue_filtered;
    double elevationAngle;
    double azimuthAngle;
};

struct SatelliteGroup {
    char name[MAX_SATELLITE_NAME_LENGTH];
    struct SatelliteData* data;
    int count;
};

void readResFiles(const char* folderPath, char* validFiles, struct SatelliteGroup* groups, int* numGroups, float interval);

#endif /* RDRESOH_H */
