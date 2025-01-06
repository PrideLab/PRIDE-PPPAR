//
// rdresoh.h
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
#ifndef RESFILE_SEARCH_H
#define RESFILE_SEARCH_H

#include <dirent.h>
#define MAX_FILENAME_LENGTH 128
#define MAX_LINE_LENGTH 1024

struct SysInfo {
    int exist;
    char sys[5];
    char freq1[5];
    char freq2[5];
};

struct FileHeader {
    char station[20];
    struct SysInfo sysInfo[10]; 
    int numSys; 
    float interval;
};

void searchResFiles(const char* folderPath, char validFiles[][MAX_FILENAME_LENGTH], int* n, struct FileHeader*);
#endif /* RESFILE_SEARCH_H */
