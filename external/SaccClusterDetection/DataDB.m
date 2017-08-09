classdef DataDB < handle%
    %  CLASS DataDB.m
    %   Reads and saves matlab variables to .mat files on disk. Simple
    %   cache implemented
    %
    % Jorge Otero-Millan, jorgeoteromillan@gmail.com 2/17/2014
    %
    %-------------------------------------------------------------------

    properties (SetAccess = private)
        folder = '';
        session = '';
    end
    
    properties (Access = private)
        cache
        USECACHE = 1;
        
        READONLY = 0;
    end
    
    methods
        
        function InitDB( this, folder, session )
            
            if ( ~ischar(folder) )
                error( 'Folder should be a string' );
            end
            if ( ~ischar(session) )
                error( 'Session should be a string' );
            end
            
            
            this.folder = folder;
            this.session = session;
            
            if ( ~exist(this.folder, 'dir') )
                error('folder does not exist');
            end
            
            if ( ~exist([this.folder '\' this.session],'dir') )
                mkdir(this.folder, this.session);
            end
        end
        
        function var = ReadVariable( this, variableName )
            var = [];
            
            try
                % if variable is in cache
                if ( isfield( this.cache, ['Session_' this.session]) && isfield( this.cache.(['Session_' this.session]), variableName ) )
                    % return variable from cache
                    var = this.cache.(['Session_' this.session]).(variableName);
                    return
                else
                    % if variable is not in cache
                    try
                        d = dir( fullfile( [this.folder '\' this.session], [variableName '.mat'] ) );
                        if isempty(d)
                            return
                        end
                        % read variable
                        dat = load(fullfile( [this.folder '\' this.session], d(1).name));
                    catch me
                        % if memory error
                        if ( isequal(me.identifier, 'MATLAB:nomem') )
                            % empty cache
                            this.cache = struct();
                            % read variable again
                            dat = load(fullfile( [this.folder '\' this.session], d(1).name));
                        else
                            rethrow(me)
                        end
                    end
                    
                    var = dat.(variableName); % TODO, change!!
                    
                    if ( this.USECACHE )
                        % add variable to cache
                        this.cache.(['Session_' this.session]).(variableName) = var;
                    end
                end
            catch me
                rethrow(me);
            end
        end
        
        function var = ReadVariableTag( this, variableName, tag )
            var = this.ReadVariable( [variableName '_' tag] );
        end
        
        function WriteVariable( this, variable, variableName )
            if ( ~this.READONLY )
                try
                    fullname = [variableName '.mat'];
                    eval([variableName ' =  variable ;']);
                    
                    save(fullfile([this.folder '\' this.session],fullname), variableName);
                    
                    if ( isfield( this.cache, ['Session_' this.session]) && isfield( this.cache.(['Session_' this.session]), variable ) )
                        this.cache.(['Session_' this.session]).(variable) = variable;
                    end
                catch me
                    rethrow(me);
                end
            else
                disp('ClusterDetection.DataDB: cannot write to the database, it is set as read only');
            end
        end

        function ClearCache( this )
            this.cache = struct();
        end
        
    end
    
    
    
end

