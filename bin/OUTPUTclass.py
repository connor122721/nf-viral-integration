#!/usr/bin/env python3

class OUTPUTfile:
    """
    Simple output file handler for tab-delimited and CSV files.
    """
    
    def __init__(self, prefix, ext, delim, columns):
        """
        Initialize output file with header.
        
        Args:
            prefix: Sample prefix for filename
            ext: File extension (e.g., '.tab', '.csv')
            delim: Delimiter character (e.g., '\t', ',')
            columns: List of column names
        """
        self.filename = prefix + ext
        self.fh = open(self.filename, "w")
        self.delimiter = delim
        
        # Write header
        self.fh.write(columns[0])
        for i in range(1, len(columns)):
            self.fh.write(delim + columns[i])
        self.fh.write('\n')
    
    def write(self, s):
        """Write string to file."""
        self.fh.write(s)
    
    def close(self):
        """Close file handle."""
        self.fh.close()