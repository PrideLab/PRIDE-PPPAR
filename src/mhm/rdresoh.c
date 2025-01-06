//
// rdresoh.c
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
// purpose   : read the header information of the residual file
//
#include "rdresoh.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


void searchResFiles(const char *folderPath, char validFiles[][MAX_FILENAME_LENGTH], int *numValidFiles, struct FileHeader* fileheader) {
    DIR *dir;
    struct dirent *entry;
    fileheader->numSys = 0;

    // Open directory
    dir = opendir(folderPath);
    if (dir == NULL) {
        perror("Error opening directory");
        return;
    }

    // Loop through files in the directory
    while ((entry = readdir(dir)) != NULL) {
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }

        // Splicing file path
        char filePath[512];
        snprintf(filePath, sizeof(filePath), "%s/%s", folderPath, entry->d_name);

        FILE *file = fopen(filePath, "r");
        if (file == NULL) {
            perror("Error opening file");
            continue;
        }

        // Read file header information
        char line[MAX_LINE_LENGTH];
        int found = 0; 
        while (fgets(line, sizeof(line), file) != NULL) {
            // Determine constellation and frequency band
            if (strstr(line, "SYS / FREQUENCY BAND") != NULL) {
                if ((*numValidFiles) == 1) {
                    struct SysInfo info;
                    info.exist = 0;
                    sscanf(line, "%s %s %s", info.sys, info.freq1, info.freq2);
                    fileheader->sysInfo[fileheader->numSys] = info;
                    (fileheader->numSys)++;
                }
            }
            else if (strstr(line, "STATION") != NULL) {
                char station[10];
                sscanf(line, "%s %*s", station);
                if ((*numValidFiles) == 0) {
                    strcpy(fileheader->station, station);
                }
                else if (strcmp(fileheader->station, station) != 0) {
                    printf("STATION fields in file header are inconsistent.");
                    exit(-1);
                }
            }

            // Determine the positioning mode of residuals
            else if (strstr(line, "POS MODE") != NULL) {
                if (strncmp(line, "Static", strlen("Static")) == 0) {
                    strcpy(validFiles[*numValidFiles], entry->d_name);
                    (*numValidFiles)++;
                    found = 1;
                } else {
                    printf("Warning: %s with wrong pos mode\n", entry->d_name);
                }

                printf("File %s format correct!\n", entry->d_name);
            }

            // Get interval
            else if (strstr(line, "INT / OBS TYPE") != NULL) {
                if ((*numValidFiles) == 1) {
                    float interval = 0.0;
                    sscanf(line, "%f %*s", &interval);
                    fileheader->interval = interval;
                }
            }
	
            else if (strstr(line, "END OF HEADER") != NULL) {
                break;
            }
        }


        fclose(file);
    }

    closedir(dir);
}
