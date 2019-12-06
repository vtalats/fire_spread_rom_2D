classdef Msgcl
  %MSGCL This class allows provides a logging utility for Matlab
  %   Allows the user to provide log level sensitive reporting for Matlab
  %   functions.  Data is output to both the log file and the console.  The
  %   log file is a running diary and is never deleted.  The format of the
  %   output is as follows:
  %
  %     hh:mm:ss |LVL | Log message string
  %
  %   where LVL is the log level for the message, hh:mm:ss is the
  %   hour:minute:second time stamp for the message, and the log message
  %   string is passed in an fprintf type syntax.  
  %
  %   When a message is passed to the class, the log level of the message is 
  %   checked against the log level of the class object.  If the the message 
  %   is less than the class then the message is output, otherwise the 
  %   message is ignored.  This allows for one to simply pass messages in
  %   the code at certain log levels and the user can control the amount of
  %   logging by simply changing the class log level
  %
  %   PUBLC
  %     Variables
  %       logLevel - This is the current log level for the class
  %
  %     Constants
  %       ALL   - -1
  %       ERR   -  2
  %       WARN  - 10
  %       PED   - 99
  %       (note these constants are used to pass the log level of a
  %       message without using specific numbers.  This gives the user the
  %       ability to group logging into broad categories)
  %
  %     Functions
  %       pmsg  - prints a formatted log message
  %               Example Usage:
  %                cls_obj.pmsg(cls_obj.ERR,'format string',<other data>);
  %       lvlck - checks a given level against the class
  %               Example Usage:
  %                 if cls_obj.lvlck(cls_obj.WARN)
  %                   plot(data);
  %                 end
  %   PRIVATE
  %     Variables
  %       filename - name of the log file
  %       fid      - file ID for the log file
  %
  %     Functions
  %       loglvlName - returns the string version of the log level number
  %
  % Author: Alan Lattimer, Virginia Tech
  % Date: April 2013
  % Versions
  %   1.0 - Original class
  %   1.1 - (April 2015) Removed logSource since it was not being used
  %
  %----------------------------------------------------------------------------
  
  properties
    logLevel = 2;
  end
  
  properties (SetAccess = private)
    ALL  = -1;
    ERR  =  2;
    WARN = 10;
    PED  = 99;
  end
  
  properties (SetAccess = private, GetAccess = private)
    fileName = '';
    fid = 1;
  end

  
  methods
    %----------------------------------------------------------------------
    % Name: Constructor
    % Description: This instantiates the class by setting the log level and
    %              opening the log file if necessary.  
    % Input: 
    %   loglevel - initial log level for the class
    %   fileName - name of the file to store the log messages.  If this is
    %              not sent then the class defaults to logging to the screen 
    %              only
    % Output:
    %   obj      - class object
    %
    function obj = Msgcl(loglevel, fileName)      
      if nargin >= 2
        obj.logLevel = loglevel;
        obj.fileName = fileName;
        obj.fid = fopen(fileName,'a');
      else
        obj.logLevel = loglevel;
      end
    end
    
    %----------------------------------------------------------------------
    % Name: pmsg
    % Description: This prints a message to the log file and to the console. 
    %              Further this determines if the message should be output
    %              by comparing its log level to the class level.  If the
    %              message <= class, then the message will be output.  The
    %              output is formatted as 
    %               hh:mm:ss |LVL | Message String
    % Input:
    %   obj      - class object
    %   msgtype  - log level for the message
    %   varargin - This is the message formatting string.  This follows the
    %              conventions for fprintf or sprintf.  See the
    %              documentation for those for more details.
    % Output:
    %   none
    %
    function pmsg(obj, msgtype, varargin)
      if obj.logLevel >= msgtype
        outstr = [datestr(now,13) ' |' obj.loglvlName(msgtype) '| ' varargin{1} '\n'];
        fprintf(outstr, varargin{2:nargin-2});
        if obj.fid ~= 1
          fprintf(obj.fid, outstr, varargin{2:nargin-2});
        end
      end
    end
    
    %----------------------------------------------------------------------
    % Name: lvlck
    % Description: This compares a log level to the class log level.  If it
    %              is less than or equal to the class then it returns true, 
    %              otherwise false.  This is useful when trying to base 
    %              other actions, such as data plot generation, on the log 
    %              level.  In this way, all conditional output can be related 
    %              to a single log level.
    % Input:
    %   obj      - class object
    %   msgtype  - log level to be checked against the class
    % Output:
    %   ck       - boolean result of the check
    %
    function ck = lvlck(obj, msgtype)
      if obj.logLevel >= msgtype
        ck = true;
      else
        ck = false;
      end
     end
    
    
    %----------------------------------------------------------------------
    % Name: Destructor
    % Description: closes the log file
    % Input:
    %   obj      - class object
    % Output:
    %   none
    %
    function delete(obj)
      fclose(obj.fid);
    end 
  
  end % Public Methods
  
  methods (Access = private)
    %----------------------------------------------------------------------
    % Name: loglvlName
    % Description: Converts special log levels to their name equivalent.
    %              If not a special log level the level is output as a left
    %              justified, space padded, 5 character string with the
    %              numeric log
    % Input:
    %   obj      - class object
    %   loglvl   - numeric log level to convert to a string
    % Output:
    %   str      - the log level number converted to a string
    %
    function str = loglvlName(obj,loglvl)
      switch loglvl
        case obj.ALL
          str = 'ALL ';
        case obj.ERR
          str = 'ERR ';
        case obj.WARN
          str = 'WARN';
        case obj.PED
          str = 'PED ';
        otherwise
          str = sprintf('%-5d',loglvl);
      end
    end

  end % Private Methods
  
end

