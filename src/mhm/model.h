//
// model.h
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
// purpose   : establishing an MHM model using filtered residual data
//
#pragma once
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include "rdresoi.h"
#include "rdresoh.h"

#define MAX_SATELLITE_GROUP 200
#define MAX_ANGLE 360
#define MAX_ELEVATION 90
#define MAX_GROUPS 26
#define MAX_BUCKETS 3000
#define ALL_FREQUENCY 18
#define MIN_POINT 1

struct Grid {
    double sumPhase, sumPseudo;
    double avgPhase, avgPseudo;
    double stdPhase, stdPseudo;
    int count, countPhase, countPseudo;
};

struct SingalFrequency {
    char singal[5];
    double frequency;
};

struct SingalFrequency* getAllFrequency() {
    struct SingalFrequency* allFrequency = malloc(ALL_FREQUENCY * sizeof(struct SingalFrequency));
    char singals[ALL_FREQUENCY][5] = {"L1", "L2", "L5", "L6", "E1", "E5a", "E6", "E5b", "E5", "B1C", "B1I", "B2a", "B3I", "B2b", "B2I", "B2", "G1", "G2"};
    double frequencys[ALL_FREQUENCY] = {1575.42, 1227.60, 1176.45, 1278.75, 1575.42, 1176.45, 1278.75, 1207.14, 1191.795, 1575.42, 1561.10, 1176.45, 1268.52, 1207.14, 1207.14, 1191.795, 1602.00, 1246.00};
    for (int i = 0; i < ALL_FREQUENCY; i++) {
        strcpy(allFrequency[i].singal, singals[i]);
        allFrequency[i].frequency = frequencys[i];
    }
    return allFrequency;    
}

void mhmHelp() {
   printf("mhm version 3.1, Wuhan University, Dec. 2024\n"); 
   printf("\n");
   printf("Usage: mhm resdir\n");
   printf("\n");
   printf("Description:\n");
   printf("mhm is  a module of PRIDE PPP-AR, is the model generator\n");
   printf("based on multipath hemispherical map\n");
   printf("\n");
   printf("Required arguments:\n");
   printf("  resdir\n");
   printf("    folder containing PRIDE PPP-AR's residual files\n");
   printf("\n");
   printf("Note: All residual files should belong to the same site\n");
   printf("\n");
   printf("Examples:\n");
   printf("  mhm resfile\n");
   printf("\n");
   printf("More details refer to PRIDE PPP-AR manual and repository\n");
   printf("  https://github.com/PrideLab/PRIDE-PPPAR/\n");
   printf("\n");
}

void getSystemName(char* sysName, char ch) {
    switch (ch) {
    case 'G':
        strcpy(sysName, "GPS");
        break;
    case 'E':
        strcpy(sysName, "GAL");
        break;
    case 'C':
        strcpy(sysName, "BDS");
        break;
    case 'R':
        strcpy(sysName, "GLO");
        break;
    case 'J':
        strcpy(sysName, "QZS");
        break;
    default:
        strcpy(sysName, "Unknown");
        break;
    }
}

char getSystemAbb(char* sysName) {
    char ch;
    if (!strcmp(sysName, "GPS")) {
        ch = 'G';
    }
    else if (!strcmp(sysName, "GAL")) {
        ch = 'E';
    }
    else if (!strcmp(sysName, "BDS")) {
        ch = 'C';
    }
    else if (!strcmp(sysName, "GLO")) {
        ch = 'R';
    }
    else if (!strcmp(sysName, "QZS")) {
        ch = 'J';
    }
    else {
        printf("Unknown system");
    }
    return ch;
}


char* extract_number_from_filename(const char* filename) {
    static char number[50]; 
    char temp[100]; 
    strcpy(temp, filename);

    char* token;
    char* rest = temp;

    while ((token = strtok(rest, "_"))) {
        if (isdigit(token[0])) {
            strcpy(number, token);
            break;
        }
        rest = NULL; 
    }

    return number;
}

void removeDuplicates(char groupCharsAll[][10], int *numGroupCharsAll) {
    int n = *numGroupCharsAll;
    int i, j, k;
    char uniqueChars[10][10];
    int uniqueCount = 0;
    for (i = 0; i < n; i++) {
        int isDuplicate = 0;
        for (j = 0; j < uniqueCount; j++) {
            if (strcmp(groupCharsAll[i], uniqueChars[j]) == 0) {
                isDuplicate = 1;
                break;
            }
        }
        if (!isDuplicate) {
            strcpy(uniqueChars[uniqueCount], groupCharsAll[i]);
            uniqueCount++;
        }
    }
    for (i = 0; i < uniqueCount; i++) {
        strcpy(groupCharsAll[i], uniqueChars[i]);
    }
    *numGroupCharsAll = uniqueCount;
}

void writeHeader(FILE* file, struct FileHeader fileheader, char validFiles[][MAX_FILENAME_LENGTH], int numValidFiles, double gridStatistic[10][2], char groupCharsAll[][10], int numGroupCharsAll) {
    fseek(file, 0, SEEK_SET);   

    fprintf(file, "MHM Model                                                   FILE TYPE\n");

    fprintf(file, "%-60sSTATION\n", fileheader.station);

    fprintf(file, "%-60dMODELING DURATION (DAY)\n", numValidFiles);

    for (int i = 0; i < numValidFiles; i++) {
        char* number = extract_number_from_filename(validFiles[i]);
        fprintf(file, "%-10s", number);
        if ((i + 1) % 5 == 0) { 
            fprintf(file, "          DAY LIST\n");
        }
    }
    if (numValidFiles % 5 != 0) { 
        for (int j = 0; j < 6 - numValidFiles % 5; j++) {
            fprintf(file, "          ");
        }
        fprintf(file, "DAY LIST\n");
    }

    char combined[11];
    for (int i = 0; i < fileheader.numSys; i++) {
        if (fileheader.sysInfo[i].exist) {
            char sysAbb;
            sysAbb = getSystemAbb(fileheader.sysInfo[i].sys);
            for (int j = 0; j < numGroupCharsAll; j++) {
                if (strchr(groupCharsAll[j], sysAbb) != NULL) {
                    snprintf(combined, sizeof(combined), "%s/%s", fileheader.sysInfo[i].freq1, fileheader.sysInfo[i].freq2);
                    fprintf(file, "%-10s%-10s%-10.4lf%-30.4lfSYS / FREQUENCY BAND / RMS (meter) / MEAN STD (meter)\n",
                        fileheader.sysInfo[i].sys, combined, gridStatistic[j][0], gridStatistic[j][1]);
                }
            }
        }
    }

    fprintf(file, "Start Field Description                                     COMMENT\n");
    fprintf(file, "AZ                                                          AZIMUTH (deg)\n");
    fprintf(file, "EL                                                          ELEVATION (deg)\n");
    fprintf(file, "P_COR                                                       PSEUDORANGE MULTIPATH CORRECTION (meter)\n");
    fprintf(file, "P_POI                                                       PSEUDORANGE MULTIPATH POINTS\n");
    fprintf(file, "P_STD                                                       PSEUDORANGE MULTIPATH STANDARD DEVIATION (meter)\n");
    fprintf(file, "P_ELI                                                       PSEUDORANGE MULTIPATH ELIMINATION RATE (%%)\n");
    fprintf(file, "L_COR                                                       CARRIER PHASE MULTIPATH CORRECTION (meter)\n");
    fprintf(file, "L_POI                                                       CARRIER PHASE MULTIPATH POINTS\n");
    fprintf(file, "L_STD                                                       CARRIER PHASE MULTIPATH STANDARD DEVIATION (meter)\n");
    fprintf(file, "L_ELI                                                       CARRIER PHASE MULTIPATH ELIMINATION RATE (%%)\n");
    fprintf(file, "End Field Description                                       COMMENT\n");
    fprintf(file, "                                                            END OF HEADER\n");
}

struct Grid** allocateGrids(int rows, int cols) {
    struct Grid** grids = (struct Grid**)malloc(rows * sizeof(struct Grid*));
    if (grids == NULL) {
        return NULL;
    }

    for (int i = 0; i < rows; i++) {
        grids[i] = (struct Grid*)malloc(cols * sizeof(struct Grid));
        if (grids[i] == NULL) {
            for (int j = 0; j < i; j++) {
                free(grids[j]);
            }
            free(grids);
            return NULL;
        }
    }

    return grids;
}

void deallocateGrids(struct Grid** grids, int rows) {
    if (grids == NULL) {
        return;
    }

    for (int i = 0; i < rows; i++) {
        free(grids[i]);
    }
    free(grids);
}

void calculateStatus(char groupCharsAll[][10], int numGroupCharsAll, struct SatelliteGroup* groups, int numGroups, struct Grid** grids, FILE* file, struct FileHeader fileheader, double gridStatistic[10][2]) {
    int rows = 360;
    int cols = 90;
    struct SatelliteData data;

    int countRecords = 0;

    int numGrid = 0;
    double gridMeanStd,gridMeanRms = 0.0;

    for (int i = 0; i < numGroupCharsAll; i++) {
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                grids[row][col].count = 0;
                grids[row][col].countPhase = 0;
                grids[row][col].countPseudo = 0;
                grids[row][col].sumPhase = 0;
                grids[row][col].sumPseudo = 0;
            }
        }

        int azimuth, elevation;
        for (int j = 0; j < numGroups; j++) {
            if (strchr(groupCharsAll[i], groups[j].name[0]) != NULL) {
                for (int k = 0; k < groups[j].count; k++) {
                    data = groups[j].data[k];
                    azimuth = (int)floor(data.azimuthAngle);
                    elevation = (int)floor(data.elevationAngle);
                    gridMeanRms = gridMeanRms + data.phaseResidue * data.phaseResidue;
                    grids[azimuth][elevation].sumPhase += data.phaseResidue;
                    grids[azimuth][elevation].sumPseudo += data.pseudoRangeResidue;
                    grids[azimuth][elevation].count += 1;
                    countRecords++;
                }
            }
        }

        for (int azimuth = 0; azimuth < rows; azimuth++) {
            for (int elevation = 0; elevation < cols; elevation++) {
                struct Grid* grid = &grids[azimuth][elevation];
                if (grid->count > 0) {
                    grid->avgPhase = grid->sumPhase / grid->count;
                    grid->avgPseudo = grid->sumPseudo / grid->count;
                }
            }
        }

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                grids[row][col].sumPhase = 0;
                grids[row][col].sumPseudo = 0;
            }
        }

        for (int j = 0; j < numGroups; j++) {
            if (strchr(groupCharsAll[i], groups[j].name[0]) != NULL) {
                for (int k = 0; k < groups[j].count; k++) {
                    data = groups[j].data[k];
                    int azimuth = (int)floor(data.azimuthAngle);
                    int elevation = (int)floor(data.elevationAngle);
                    grids[azimuth][elevation].sumPhase += pow(data.phaseResidue - grids[azimuth][elevation].avgPhase, 2);
                    grids[azimuth][elevation].sumPseudo += pow(data.pseudoRangeResidue - grids[azimuth][elevation].avgPseudo, 2);
                }
            }
        }

        for (int azimuth = 0; azimuth < rows; azimuth++) {
            for (int elevation = 0; elevation < cols; elevation++) {
                struct Grid* grid = &grids[azimuth][elevation];
                if (grid->count > 0) {
                    grid->stdPhase = sqrt(grid->sumPhase / grid->count);
                    grid->stdPseudo = sqrt(grid->sumPseudo / grid->count);
                }
            }
        }

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                grids[row][col].sumPhase = 0;
                grids[row][col].sumPseudo = 0;
            }
        }


        for (int j = 0; j < numGroups; j++) {
            if (strchr(groupCharsAll[i], groups[j].name[0]) != NULL) {
                for (int k = 0; k < groups[j].count; k++) {

                    data = groups[j].data[k];
                    int azimuth = (int)floor(data.azimuthAngle);
                    int elevation = (int)floor(data.elevationAngle);

                    if (fabs(grids[azimuth][elevation].avgPhase - data.phaseResidue) <= 2 * grids[azimuth][elevation].stdPhase) {
                        grids[azimuth][elevation].sumPhase += data.phaseResidue;
                        grids[azimuth][elevation].countPhase += 1;
                    }
                    if (fabs(grids[azimuth][elevation].avgPseudo - data.pseudoRangeResidue) <= 2 * grids[azimuth][elevation].stdPseudo) {
                        grids[azimuth][elevation].sumPseudo += data.pseudoRangeResidue;
                        grids[azimuth][elevation].countPseudo += 1;
                    }

                }
            }
        }

        for (int azimuth = 0; azimuth < rows; azimuth++) {
            for (int elevation = 0; elevation < cols; elevation++) {
                struct Grid* grid = &grids[azimuth][elevation];
                if (grid->count > 0) {
                    grid->avgPhase = grid->sumPhase / grid->countPhase;
                    grid->avgPseudo = grid->sumPseudo / grid->countPseudo;
                }
            }
        }

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                grids[row][col].sumPhase = 0;
                grids[row][col].sumPseudo = 0;
            }
        }

        for (int j = 0; j < numGroups; j++) {
            if (strchr(groupCharsAll[i], groups[j].name[0]) != NULL) {
                for (int k = 0; k < groups[j].count; k++) {
                    data = groups[j].data[k];
                    int azimuth = (int)floor(data.azimuthAngle);
                    int elevation = (int)floor(data.elevationAngle);

                    if (fabs(grids[azimuth][elevation].avgPhase - data.phaseResidue) <= 2 * grids[azimuth][elevation].stdPhase) {
                        grids[azimuth][elevation].sumPhase += pow(data.phaseResidue - grids[azimuth][elevation].avgPhase, 2);
                    }

                    if (fabs(grids[azimuth][elevation].avgPseudo - data.pseudoRangeResidue) <= 2 * grids[azimuth][elevation].stdPseudo) {
                        grids[azimuth][elevation].sumPseudo += pow(data.pseudoRangeResidue - grids[azimuth][elevation].avgPseudo, 2);
                    }
                }
            }
        }

        for (int azimuth = 0; azimuth < rows; azimuth++) {
            for (int elevation = 0; elevation < cols; elevation++) {
                struct Grid* grid = &grids[azimuth][elevation];
                if (grid->count > 0) {
                    grid->stdPhase = sqrt(grid->sumPhase / grid->countPhase);
                    grid->stdPseudo = sqrt(grid->sumPseudo / grid->countPseudo);
                    gridMeanStd += grid->stdPhase; 
                    numGrid++;
                }
            }
        }
        gridMeanStd = gridMeanStd / numGrid;
        gridMeanRms = gridMeanRms / countRecords;
        gridMeanRms = sqrt(gridMeanRms);
        gridStatistic[i][0] = gridMeanRms;
        gridStatistic[i][1] = gridMeanStd;     

        int length = strlen(groupCharsAll[i]);

        fprintf(file, "                                                            START OF ");
        for (int j= 0; j < length; j++) {
            char sysName[20]; 
            if (groupCharsAll[i][j]== 'G') {
                strcpy(sysName, "GPS");
            }
            else if (groupCharsAll[i][j] == 'E') {
                strcpy(sysName, "GAL");
            }
            else if (groupCharsAll[i][j] == 'C') {
                strcpy(sysName, "BDS");
            }
            else if (groupCharsAll[i][j] == 'R') {
                strcpy(sysName, "GLO");
            }
            else if (groupCharsAll[i][j] == 'J') {
                strcpy(sysName, "QZS");
            }
            else {
                strcpy(sysName, "Unknown");
            }

            int fileheaderIndex = 0;
            for (int k = 0; k < fileheader.numSys; k++) {
                if (!strcmp(sysName, fileheader.sysInfo[k].sys)) {
                    fileheaderIndex = k;
                }
            }

            fprintf(file, "%s %s/%s ",
                fileheader.sysInfo[fileheaderIndex].sys, fileheader.sysInfo[fileheaderIndex].freq1, fileheader.sysInfo[fileheaderIndex].freq2);
        }
        fprintf(file,"\n");
        fprintf(file, "   AZ   EL     P_COR  P_POI     P_STD     P_ELI     L_COR  L_POI     L_STD     L_ELI\n");
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                struct Grid* grid = &grids[row][col];
                if (grid->count > MIN_POINT) {
                    fprintf(file, "%5d", row);
                    fprintf(file, "%5d", col);

                    fprintf(file, "%10.4lf", grid->avgPseudo);
                    fprintf(file, "%7d", grid->countPseudo);
                    fprintf(file, "%10.4lf", grid->stdPseudo);
                    fprintf(file, "%10.4lf", (1 - grid->countPseudo * 1.0 / (grid->count * 1.0)) * 100);

                    fprintf(file, "%10.4lf", grid->avgPhase);
                    fprintf(file, "%7d", grid->countPhase);
                    fprintf(file, "%10.4lf", grid->stdPhase);
                    fprintf(file, "%10.4lf", (1 - grid->countPhase * 1.0 / (grid->count * 1.0)) * 100);
                    fprintf(file, "\n"); 

                }
            }
        }
        fprintf(file, "                                                            END OF ");
        for (int j= 0; j < length; j++) {
            char sysName[20];
            if (groupCharsAll[i][j] == 'G') {
                strcpy(sysName, "GPS");
            }
            else if (groupCharsAll[i][j] == 'E') {
                strcpy(sysName, "GAL");
            }
            else if (groupCharsAll[i][j] == 'C') {
                strcpy(sysName, "BDS");
            }
            else if (groupCharsAll[i][j] == 'R') {
                strcpy(sysName, "GLO");
            }
            else if (groupCharsAll[i][j] == 'J') {
                strcpy(sysName, "QZS");
            }
            else {
                strcpy(sysName, "Unknown");
            }

            int fileheaderIndex = 0;
            for (int k = 0; k < fileheader.numSys; k++) {
                if (!strcmp(sysName, fileheader.sysInfo[k].sys)) {
                    fileheaderIndex = k;
                }
            }

            fprintf(file, "%s %s/%s ",
                fileheader.sysInfo[fileheaderIndex].sys, fileheader.sysInfo[fileheaderIndex].freq1, fileheader.sysInfo[fileheaderIndex].freq2);
        }
        fprintf(file,"\n");

    }

// printf("countAllRecords: %d\n", countRecords);
}
