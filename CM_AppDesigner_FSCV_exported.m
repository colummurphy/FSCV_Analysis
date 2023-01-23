classdef CM_AppDesigner_FSCV_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure      matlab.ui.Figure
        Menu          matlab.ui.container.Menu
        LoadMenu      matlab.ui.container.Menu
        ColorCheck4   matlab.ui.control.CheckBox
        ColorCheck3   matlab.ui.control.CheckBox
        ColorCheck2   matlab.ui.control.CheckBox
        ColorCheck1   matlab.ui.control.CheckBox
        DataTable     matlab.ui.control.Table
        TitleText     matlab.ui.control.EditField
        CvPlot        matlab.ui.control.UIAxes
        DopaminePlot  matlab.ui.control.UIAxes
        ColorPlot4    matlab.ui.control.UIAxes
        ColorPlot3    matlab.ui.control.UIAxes
        ColorPlot2    matlab.ui.control.UIAxes
        ColorPlot1    matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        % Property % Description
        % Property2 % Description


    end
    
    methods (Access = private)

        function startScript(app)

            % make app invisible until uigetfile is displayed
            app.UIFigure.Visible = 'off'; 

            % get the file here
            [selectedFileName,selectedPathName] = ... 
                        uigetfile('*.*', 'Select file');

            app.UIFigure.Visible = 'on'; 
            
            % create color plot array
            appDesColorPlots = [app.ColorPlot1, app.ColorPlot2, ... 
                                    app.ColorPlot3, app.ColorPlot4];
            
            % create check box array
            checkBoxes = [app.ColorCheck1, app.ColorCheck2,...
                          app.ColorCheck3, app.ColorCheck4];
            
            checkNums = [1, 2, 3, 4];

            % pass the app designer objects to the main fscv function
            CM_plotEphysFSCVmulti(app.UIFigure, selectedFileName,... 
                   selectedPathName, appDesColorPlots, app.DopaminePlot,... 
                   app.TitleText, checkBoxes, checkNums);

            % align the dopamine plot to the first color plot
            app.DopaminePlot.InnerPosition(1) = app.ColorPlot1.InnerPosition(1);
            app.DopaminePlot.InnerPosition(3) = app.ColorPlot1.InnerPosition(3);
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % call the start script
            startScript(app);
        end

        % Window key press function: UIFigure
        function UIFigureWindowKeyPress(app, event)
            keyPressed = event.Key;

            appDesColorPlots = [app.ColorPlot1, app.ColorPlot2, ... 
                                    app.ColorPlot3, app.ColorPlot4];
            
            CM_buttonpress(keyPressed, appDesColorPlots, ...
                       app.DopaminePlot, app.TitleText);
        end

        % Window button down function: UIFigure
        function UIFigureWindowButtonDown(app, event)
           
            % figure -> axes -> image
            figureObj = event.Source;
            childObj = event.Source.CurrentObject;
            
            % Values - 'normal' (left), 'alt' (right), 'extend' (center)
            typeOfClick = figureObj.SelectionType;

            % [x, y] pixels
            clickPosOnFigure = figureObj.CurrentPoint;

            % if there is a child object and it is an image
            if ~isempty(childObj) && strcmp(childObj.Type, 'image')
           
                % get the plot position from the parent (axis object)
                plotPosOnFigure = childObj.Parent.InnerPosition;

                % disp width + height
                % disp(["width of plot: " plotPosOnFigure(3)]);
                % disp(["height of plot: " plotPosOnFigure(4)]);

                % [x, y] - [left bottom width height]       
                clickPosOnPlot = ...
                 [clickPosOnFigure(1) - plotPosOnFigure(1),...
                  clickPosOnFigure(2) - plotPosOnFigure(2)];
                
                appDesColorPlots = [app.ColorPlot1, app.ColorPlot2, ... 
                                    app.ColorPlot3, app.ColorPlot4];

                % call mouse handler function
                CM_mouseclick(typeOfClick, plotPosOnFigure, clickPosOnPlot, ...
                   appDesColorPlots, app.CvPlot, app.DataTable, app.DopaminePlot);

            end
        end

        % Value changed function: ColorCheck1
        function ColorCheck1ValueChanged(app, event)
            
            % Color Plot 1  
            % value from checkbox (0 or 1)
            checkNum = 1;
            checkVal = app.ColorCheck1.Value;            
            
            appDesColorPlots = [app.ColorPlot1, app.ColorPlot2, ... 
                              app.ColorPlot3, app.ColorPlot4];
            
            CM_checkbox(checkNum, checkVal, appDesColorPlots, ...
                        app.DopaminePlot);
            
        end

        % Value changed function: ColorCheck2
        function ColorCheck2ValueChanged(app, event)
            
            % Color Plot 1
            % value from checkbox (0 or 1)
            checkNum = 2;
            checkVal = app.ColorCheck2.Value;            
            
            
            appDesColorPlots = [app.ColorPlot1, app.ColorPlot2, ... 
                              app.ColorPlot3, app.ColorPlot4];
            
            CM_checkbox(checkNum, checkVal, appDesColorPlots, ...
                        app.DopaminePlot);

        end

        % Value changed function: ColorCheck3
        function ColorCheck3ValueChanged(app, event)
            
            % Color Plot 3
            % value from checkbox (0 or 1)
            checkNum = 3;
            checkVal = app.ColorCheck3.Value;            
            
            appDesColorPlots = [app.ColorPlot1, app.ColorPlot2, ... 
                              app.ColorPlot3, app.ColorPlot4];
            
            CM_checkbox(checkNum, checkVal, appDesColorPlots, ...
                        app.DopaminePlot);

        end

        % Value changed function: ColorCheck4
        function ColorCheck4ValueChanged(app, event)

            % Color Plot 4  
            % value from checkbox (0 or 1)
            checkNum = 4;
            checkVal = app.ColorCheck4.Value;            
            
            appDesColorPlots = [app.ColorPlot1, app.ColorPlot2, ... 
                              app.ColorPlot3, app.ColorPlot4];
            
            CM_checkbox(checkNum, checkVal, appDesColorPlots, ...
                        app.DopaminePlot);

        end

        % Menu selected function: LoadMenu
        function LoadMenuSelected(app, event)
           
            appDesColorPlots = [app.ColorPlot1, app.ColorPlot2, ... 
                                    app.ColorPlot3, app.ColorPlot4];
            
            CM_menuselect(appDesColorPlots, app.DopaminePlot, app.TitleText);
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [20 20 1250 700];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.WindowButtonDownFcn = createCallbackFcn(app, @UIFigureWindowButtonDown, true);
            app.UIFigure.WindowKeyPressFcn = createCallbackFcn(app, @UIFigureWindowKeyPress, true);

            % Create Menu
            app.Menu = uimenu(app.UIFigure);
            app.Menu.Text = 'Menu';

            % Create LoadMenu
            app.LoadMenu = uimenu(app.Menu);
            app.LoadMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadMenuSelected, true);
            app.LoadMenu.Text = 'Load';

            % Create ColorPlot1
            app.ColorPlot1 = uiaxes(app.UIFigure);
            ylabel(app.ColorPlot1, 'Voltage (V)')
            zlabel(app.ColorPlot1, 'Z')
            app.ColorPlot1.XLimitMethod = 'tight';
            app.ColorPlot1.YLimitMethod = 'tight';
            app.ColorPlot1.XTick = [];
            app.ColorPlot1.XTickLabel = '';
            app.ColorPlot1.YTick = [];
            app.ColorPlot1.Tag = 'ColorPlot1';
            app.ColorPlot1.Position = [90 539 750 128];

            % Create ColorPlot2
            app.ColorPlot2 = uiaxes(app.UIFigure);
            ylabel(app.ColorPlot2, 'Voltage (V)')
            zlabel(app.ColorPlot2, 'Z')
            app.ColorPlot2.XLimitMethod = 'tight';
            app.ColorPlot2.YLimitMethod = 'tight';
            app.ColorPlot2.XTick = [];
            app.ColorPlot2.XTickLabel = '';
            app.ColorPlot2.YTick = [];
            app.ColorPlot2.Tag = 'ColorPlot2';
            app.ColorPlot2.Position = [90 406 750 128];

            % Create ColorPlot3
            app.ColorPlot3 = uiaxes(app.UIFigure);
            ylabel(app.ColorPlot3, 'Voltage (V)')
            zlabel(app.ColorPlot3, 'Z')
            app.ColorPlot3.XLimitMethod = 'tight';
            app.ColorPlot3.YLimitMethod = 'tight';
            app.ColorPlot3.XTick = [];
            app.ColorPlot3.XTickLabel = '';
            app.ColorPlot3.YTick = [];
            app.ColorPlot3.Tag = 'ColorPlot3';
            app.ColorPlot3.Position = [90 273 750 128];

            % Create ColorPlot4
            app.ColorPlot4 = uiaxes(app.UIFigure);
            ylabel(app.ColorPlot4, 'Voltage (V)')
            zlabel(app.ColorPlot4, 'Z')
            app.ColorPlot4.XLimitMethod = 'tight';
            app.ColorPlot4.YLimitMethod = 'tight';
            app.ColorPlot4.XTick = [];
            app.ColorPlot4.XTickLabel = '';
            app.ColorPlot4.YTick = [];
            app.ColorPlot4.Tag = 'ColorPlot4';
            app.ColorPlot4.Position = [90 140 750 128];

            % Create DopaminePlot
            app.DopaminePlot = uiaxes(app.UIFigure);
            ylabel(app.DopaminePlot, 'DM')
            zlabel(app.DopaminePlot, 'Z')
            app.DopaminePlot.XTick = [];
            app.DopaminePlot.XTickLabel = '';
            app.DopaminePlot.TickDir = 'out';
            app.DopaminePlot.Position = [90 5 750 130];

            % Create CvPlot
            app.CvPlot = uiaxes(app.UIFigure);
            title(app.CvPlot, 'CV at selected timepoint')
            xlabel(app.CvPlot, 'E_a (V)')
            ylabel(app.CvPlot, 'I (nA)')
            zlabel(app.CvPlot, 'Z')
            app.CvPlot.XLim = [-0.4 1.3];
            app.CvPlot.XTick = [-0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.3];
            app.CvPlot.XTickLabel = {'-0.4'; '-0.2'; '0'; '0.2'; '0.4'; '0.6'; '0.8'; '1'; '1.3'};
            app.CvPlot.Position = [855 57 384 280];

            % Create TitleText
            app.TitleText = uieditfield(app.UIFigure, 'text');
            app.TitleText.Position = [855 630 376 22];

            % Create DataTable
            app.DataTable = uitable(app.UIFigure);
            app.DataTable.ColumnName = {'Ch1'; 'Ch2'; 'Ch3'; 'Ch4'};
            app.DataTable.RowName = {};
            app.DataTable.FontSize = 10;
            app.DataTable.Position = [855 355 376 244];

            % Create ColorCheck1
            app.ColorCheck1 = uicheckbox(app.UIFigure);
            app.ColorCheck1.ValueChangedFcn = createCallbackFcn(app, @ColorCheck1ValueChanged, true);
            app.ColorCheck1.Text = 'Plot1';
            app.ColorCheck1.Position = [20 584 50 22];

            % Create ColorCheck2
            app.ColorCheck2 = uicheckbox(app.UIFigure);
            app.ColorCheck2.ValueChangedFcn = createCallbackFcn(app, @ColorCheck2ValueChanged, true);
            app.ColorCheck2.Text = 'Plot2';
            app.ColorCheck2.Position = [20 454 50 22];

            % Create ColorCheck3
            app.ColorCheck3 = uicheckbox(app.UIFigure);
            app.ColorCheck3.ValueChangedFcn = createCallbackFcn(app, @ColorCheck3ValueChanged, true);
            app.ColorCheck3.Text = 'Plot3';
            app.ColorCheck3.Position = [20 325 50 22];

            % Create ColorCheck4
            app.ColorCheck4 = uicheckbox(app.UIFigure);
            app.ColorCheck4.ValueChangedFcn = createCallbackFcn(app, @ColorCheck4ValueChanged, true);
            app.ColorCheck4.Text = 'Plot4';
            app.ColorCheck4.Position = [20 194 50 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = CM_AppDesigner_FSCV_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end