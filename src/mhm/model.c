//
// model.c
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
// purpose   : establishing an MHM model using residual data
//
#include <stdio.h>
#include "model.h"

int main(int argc, char *argv[]) {
    struct FileHeader fileheader; //File header information
    char validFiles[MAX_FILENAME_LENGTH][MAX_FILENAME_LENGTH]; // Valid file name
    int numValidFiles = 0; // Number of valid files
    
    if(argc != 2) {
        mhmHelp();
        exit(1);
    }
    
    // Search resfile folder
    searchResFiles(argv[1], validFiles, &numValidFiles, &fileheader);

    printf("fileheader.station:%s\n", fileheader.station);
    printf("fileheader.interval:%f\n", fileheader.interval);

    char outputFileName[24];
    snprintf(outputFileName, sizeof(outputFileName), "mhm_%s", fileheader.station);
    FILE* file; 
    file = fopen(outputFileName, "w");
    if (file == NULL) {
        printf("Unable to open file.\n");
        return -1;
    }

    struct SatelliteGroup* allGroups = (struct SatelliteGroup*)malloc(MAX_SATELLITE_GROUP * 10 * sizeof(struct SatelliteGroup));
    if (allGroups == NULL) {
        fprintf(stderr, "Failed to allocate memory for allGroups\n");
        return -1; 
    }
    int allGroupsCount = 0;

    int rows = 360;
    int cols = 90;
    struct Grid** grids = allocateGrids(rows, cols);
    if (grids == NULL) {
        printf("Failed to allocate memory for allGroups\n");
        return -1;
    }

    //Loop through each file
    for (int fileInd = 0; fileInd < numValidFiles; fileInd++) {
        printf("%s\n", validFiles[fileInd]);
        struct SatelliteGroup groups[MAX_SATELLITE_GROUP] = { 0 };
        int numGroups = 0;
        //Group the returned data by satellite
        readResFiles(argv[1], validFiles[fileInd], groups, &numGroups, fileheader.interval);
        printf("Num groups:%d\n\n", numGroups);

        for (int i = 0; i < numGroups; i++) {
            memcpy(&allGroups[allGroupsCount], &groups[i], sizeof(struct SatelliteGroup));
            allGroupsCount++;
        }

    }

    printf("count all groups: %d\n", allGroupsCount);
    char groupChars[MAX_GROUPS] = { 0 };
    char groupCharsAll[MAX_GROUPS][10];
    memset(groupCharsAll, 0, sizeof(groupCharsAll));
    int numGroupChars = 0;
    struct SingalFrequency *allFrequency = getAllFrequency();
    for (int i = 0; i < allGroupsCount; i++) {
        if (strchr(groupChars, allGroups[i].name[0]) == NULL) {
            groupChars[numGroupChars] = allGroups[i].name[0];
            char sysName[20];
            getSystemName(sysName, allGroups[i].name[0]);
            for (int j = 0; j < fileheader.numSys; j++) {
                if (!strcmp(sysName, fileheader.sysInfo[j].sys)) {
                    fileheader.sysInfo[j].exist = 1;
                }
            }
            numGroupChars++;
        }
    }

    for (int i = 0; i < numGroupChars; i++) {
        char sysName[20];
        int  sysIndex;
        getSystemName(sysName, groupChars[i]);
        for (int j = 0; j < fileheader.numSys; j++) {
            if (!strcmp(sysName, fileheader.sysInfo[j].sys)) {
                sysIndex = j;
            }
        }
        for (int j = 0; j < fileheader.numSys; j++) {
            double freqIndex1 = 0.0;
            double freqIndex2 = 0.0;
            double freqCompare1 = 1.0;
            double freqCompare2 = 1.0;
            for (int k = 0; k < ALL_FREQUENCY; k++) {
                if (!strcmp(fileheader.sysInfo[sysIndex].freq1, allFrequency[k].singal)) {
                    freqIndex1 = allFrequency[k].frequency;
                }
                if (!strcmp(fileheader.sysInfo[sysIndex].freq2, allFrequency[k].singal)) {
                    freqIndex2 = allFrequency[k].frequency;
                }
                if (!strcmp(fileheader.sysInfo[j].freq1, allFrequency[k].singal)) {
                    freqCompare1 = allFrequency[k].frequency;
                }
                if (!strcmp(fileheader.sysInfo[j].freq2, allFrequency[k].singal)) {
                    freqCompare2 = allFrequency[k].frequency;
                }
            }
            if (freqIndex1 == freqCompare1 && freqIndex2 == freqCompare2) { 
                char sysAbb;
                sysAbb = getSystemAbb(fileheader.sysInfo[j].sys);
                int len = strlen(groupCharsAll[i]);
                if (len + 1 < 10) {
                    groupCharsAll[i][len] = sysAbb;
                    groupCharsAll[i][len + 1] = '\0';
                } 
            }
        }
    }

    free(allFrequency);

    int numGroupCharsAll = numGroupChars;
    removeDuplicates(groupCharsAll, &numGroupCharsAll);

    double gridStatistic[10][2];
   
    calculateStatus(groupCharsAll, numGroupCharsAll, allGroups, allGroupsCount, grids, file, fileheader, gridStatistic); 
    writeHeader(file, fileheader, validFiles, numValidFiles, gridStatistic, groupCharsAll, numGroupCharsAll);
    calculateStatus(groupCharsAll, numGroupCharsAll, allGroups, allGroupsCount, grids, file, fileheader, gridStatistic);

    printf("Program execution completed.\n");
    fclose(file);
    deallocateGrids(grids, rows);
    return 0;
}

