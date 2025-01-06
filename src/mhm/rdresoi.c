//
// rdresoi.c
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
#include "rdresoi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Read all satellite data files in the resfile folder
void readResFiles(const char* folderPath, char* validFiles, struct SatelliteGroup* groups, int* numGroups, float interval) {
    char filePath[MAX_FILENAME_LENGTH];
    snprintf(filePath, sizeof(filePath), "%s/%s", folderPath, validFiles);

    FILE* file = fopen(filePath, "r");
    if (file == NULL) {
        perror("Error opening file");
    }

    char line[MAX_LINE_LENGTH];
    int startReadingData = 0;
    float timeInSeconds;
    int numRecords = 0;

    while (fgets(line, sizeof(line), file) != NULL) {
        if (strncmp(line, "TIM", 3) == 0) {
            sscanf(line, "%*s %*s %*s %*s %*s %*s %*s %*s %f", &timeInSeconds);
            startReadingData = 1;
        }
        else if (startReadingData) {
            numRecords += 1;
            // Analyze satellite data
            char satelliteName[MAX_SATELLITE_NAME_LENGTH];
            struct SatelliteData data;
            sscanf(line, "%s %lf %lf %*s %*s %*s %lf %lf",
                satelliteName, &data.phaseResidue,
                &data.pseudoRangeResidue, &data.elevationAngle,
                &data.azimuthAngle);
            // Naturalization azimuth
            if (data.azimuthAngle < 0) {
                data.azimuthAngle += 360;
            }

            data.timeInSeconds = timeInSeconds;

            // Record satellite name
            int isNewSatellite = 1;
            for (int j = 0; j < *numGroups; ++j) {
                if (strcmp(groups[j].name, satelliteName) == 0) {
                    groups[j].data[groups[j].count++] = data;
                    isNewSatellite = 0;
                    break;
                }
            }

            if (isNewSatellite && (*numGroups) < MAX_SATELLITES) {
                strcpy(groups[*numGroups].name, satelliteName);
				int numEpoch = (int)(86400 / interval);
                groups[*numGroups].data = malloc(sizeof(struct SatelliteData) * numEpoch); 
                if (groups[*numGroups].data == NULL) {
                    fprintf(stderr, "Failed to allocate memory for groups[%d].data\n", *numGroups);
                    exit(EXIT_FAILURE); 
                }
                groups[*numGroups].count = 1;
                groups[*numGroups].data[0] = data;
                (*numGroups)++;
            }
        }
    }

    printf("numRecords: %d\n", numRecords);

    fclose(file);

}
