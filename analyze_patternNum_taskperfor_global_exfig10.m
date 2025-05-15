clear;clc
addpath(genpath('functions'));
hemisphere = {'lh','rh'};
parcellation = 'HCPex_360';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'

[mask_L,label_L,oi_x_L,oi_y_L,V_mask_L] = masklabel(parcellation,'L');
[mask_R,label_R,oi_x_R,oi_y_R,V_mask_R] = masklabel(parcellation,'R');

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R

Condition = {'WM', 'RELATIONAL'};
PE = {'RL','LR'};


figure('color','w')
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,18,22.5]);

t0 = tiledlayout(4,1,"TileSpacing","compact","Padding",'compact');
for cdn = 1:2
    for pe = 1:2

        switch Condition{cdn}
            case 'WM'
                TaskName = 'Working Memory';
            case 'RELATIONAL'
                TaskName = 'Relational Processing';
        end


        switch PE{pe}
            case 'RL'
                Run = '1';
            case 'LR'
                Run = '2';
        end
        load(['data/results/info_task/' Condition{cdn} '_' PE{pe} '_info.mat']);
        load(['data/results/' Condition{cdn} '_' PE{pe} '_Sublist.mat']);

        % global cortex level
        if length(hemisphere)==2
            Num_L = load(['data/results/numPattern_results_new/global/' hemisphere{1} '_numPattern_individual_' Condition{cdn} '_' PE{pe} '.mat']);
            Num_R = load(['data/results/numPattern_results_new/global/' hemisphere{2} '_numPattern_individual_' Condition{cdn} '_' PE{pe} '.mat']);
            NumSource = Num_L.NumSource + Num_R.NumSource;
            NumSink = Num_L.NumSink + Num_R.NumSink;
            NumCCVortex = Num_L.NumCCVortex + Num_R.NumCCVortex;
            NumCVortex = Num_L.NumCVortex + Num_R.NumCVortex;
        else
            load(['data/results/numPattern_results_new/global/' hemisphere{1} '_numPattern_individual_' Condition{cdn} '_' PE{pe} '.mat']);
        end

        if strcmp(Condition{cdn},'WM')
            Acc(Acc==0)=nan;
        end
        Subind=find(~isnan(Acc(:,6)));% 1:size(NumSource,1);
        T = 5:size(NumSource,2)-4;
        T(T>size(NumSource,2))=[];
        NumSource = NumSource(Subind,T);
        NumSink = NumSink(Subind,T);
        NumCCVortex = NumCCVortex(Subind,T);
        NumCVortex = NumCVortex(Subind,T);

        Acc = Acc(Subind,~isnan(Acc(1,:)));
        [transformed_data, lambda] = boxcox(mean(Acc,2,'omitnan'));
        % Acc = log(Acc ./ (1 - Acc));
        RT = RT(Subind,~isnan(RT(1,:)));
        Threshold = 0;

        SubjectList = SubjectList(Subind);

        behav_info = readtable('data/fMRI/unrestricted_qchyyyyyy_11_22_2021_0_28_6.csv');
        behav_info = behav_info(ismember(behav_info.Subject,str2num(cell2mat(SubjectList))),:);

        sub_gender = behav_info.Gender;
        sub_gender_bi = zeros(size(sub_gender));
        sub_gender_bi(cellfun(@(x) x=='M',sub_gender))=1;
        sub_age = behav_info.Age;
        sub_age = cellfun(@(x) split(x,'-'), sub_age,'UniformOutput',false);
        [sub_age{cell2mat(cellfun(@(x) length(x)<2,sub_age,'UniformOutput',false))}] = deal({'36';'36'});
        sub_age = cell2mat(cellfun(@(x) (str2double(x{1})+str2double(x{2}))/2, sub_age,'UniformOutput',false));

        pattern_names_old = {'Source','Sink','CCVortex','CVortex'};
        pattern_names = {'Diverging','Converging','CC-spiral','C-spiral'};

        nexttile(t0)
        set(gca,'Visible','off');

        t1 = tiledlayout(t0,2,5,"TileSpacing","compact","Padding",'compact');
        t1.Layout.Tile = 2*(cdn-1)+pe;
        
        if 2*(cdn-1)+pe == 4; xlabel(t1,'The number of Patterns'); end
        % title(t1,[TaskName ' Run ' Run]);
        set(gca,'FontName','Arial','FontUnits','points','FontSize',6);
        %
        YYLim = [mean(transformed_data)-3*std(transformed_data) mean(transformed_data)+3*std(transformed_data)];
        t = []; p = [];
        for i_pattern = 1:4
            nexttile(t1)
            eval(['x=mean(Num' pattern_names_old{i_pattern} ',2);'])
            Y=transformed_data;
            X_design = [x sub_age sub_gender_bi];

            [b,stats] = robustfit(X_design,Y);
            t(i_pattern,1) = stats.t(2);
            P(i_pattern,1) = stats.p(2);

            XXLim = [mean(x)-3*std(x), mean(x)+3*std(x)];
            x_fit = linspace(XXLim(1), XXLim(2), 100)';
            X_pred = [x_fit, mean(sub_age)*ones(size(x_fit)), mean(sub_gender_bi)*ones(size(x_fit))];
            y_fit = [ones(size(x_fit)), X_pred] * b;

            h3=scatter(x,Y,[],[0.6 0.6 0.6],"o",'filled','LineWidth',0.05);alpha(h3,0.5);
            hold on
            plot(x_fit,y_fit,'b','linewidth',1);
            xlim(XXLim);ylim(YYLim);
            title(pattern_names{i_pattern});
            if i_pattern==1;ylabel('Accuracy (boxcox)');end

            text(0.95*XXLim(1)+0.05*XXLim(2),0.95*YYLim(2)+0.05*YYLim(1),{['P = ' num2str(stats.p(2),'%.2e')]});

            set(gca,'FontName','Arial','FontUnits','points','FontSize',6);
        end
        
        nexttile(t1)
        x=mean(NumSource,2)+mean(NumSink,2)+mean(NumCCVortex,2)+mean(NumCVortex,2);
        Y=transformed_data;
        X_design = [x sub_age sub_gender_bi];

        [b,stats] = robustfit(X_design,Y);
        t(5,1) = stats.t(2);
        P(5,1) = stats.p(2);

        XXLim = [mean(x)-3*std(x), mean(x)+3*std(x)];
        x_fit = linspace(XXLim(1), XXLim(2), 100)';
        X_pred = [x_fit, mean(sub_age)*ones(size(x_fit)), mean(sub_gender_bi)*ones(size(x_fit))];
        y_fit = [ones(size(x_fit)), X_pred] * b;

        h3=scatter(x,Y,[],[0.6 0.6 0.6],"o",'filled','LineWidth',0.05);alpha(h3,0.5);
        hold on
        plot(x_fit,y_fit,'b','linewidth',1);
        xlim(XXLim);ylim(YYLim);
        title('All');
        text(0.95*XXLim(1)+0.05*XXLim(2),0.95*YYLim(2)+0.05*YYLim(1),{['P = ' num2str(stats.p(2),'%.2e')]});
        set(gca,'FontName','Arial','FontUnits','points','FontSize',6);

        for i_pattern = 1:4
            nexttile(t1)
            
            eval(['x=mean(Num' pattern_names_old{i_pattern} ',2);'])
            Y=mean(RT,2,'omitnan');
            X_design = [x sub_age sub_gender_bi];

            [b,stats] = robustfit(X_design,Y);
            t(i_pattern,2) = stats.t(2);
            P(i_pattern,2) = stats.p(2);

            XXLim = [mean(x)-3*std(x), mean(x)+3*std(x)];
            YYLim = [mean(Y)-3*std(Y) mean(Y)+3*std(Y)];

            x_fit = linspace(XXLim(1), XXLim(2), 100)';
            X_pred = [x_fit, mean(sub_age)*ones(size(x_fit)), mean(sub_gender_bi)*ones(size(x_fit))];
            y_fit = [ones(size(x_fit)), X_pred] * b;

            h3=scatter(x,Y,[],[0.6 0.6 0.6],"o",'filled','LineWidth',0.05);alpha(h3,0.5);
            hold on
            plot(x_fit,y_fit,'b','linewidth',1);
            xlim(XXLim);ylim(YYLim);
            title(pattern_names{i_pattern});
            if i_pattern==1;ylabel('Reaction time');end

            text(0.95*XXLim(1)+0.05*XXLim(2),0.95*YYLim(2)+0.05*YYLim(1),{['P = ' num2str(stats.p(2),'%.2e')]});
            set(gca,'FontName','Arial','FontUnits','points','FontSize',6);
        end

        nexttile(t1)
        
        x=mean(NumSource,2)+mean(NumSink,2)+mean(NumCCVortex,2)+mean(NumCVortex,2);
        Y=mean(RT,2,'omitnan');
        X_design = [x sub_age sub_gender_bi];

        [b,stats] = robustfit(X_design,Y);
        t(5,2) = stats.t(2);
        P(5,2) = stats.p(2);

        XXLim = [mean(x)-3*std(x), mean(x)+3*std(x)];
        YYLim = [mean(Y)-3*std(Y) mean(Y)+3*std(Y)];

        x_fit = linspace(XXLim(1), XXLim(2), 100)';
        X_pred = [x_fit, mean(sub_age)*ones(size(x_fit)), mean(sub_gender_bi)*ones(size(x_fit))];
        y_fit = [ones(size(x_fit)), X_pred] * b;

        h3=scatter(x,Y,[],[0.6 0.6 0.6],"o",'filled','LineWidth',0.05);alpha(h3,0.5);
        hold on
        h2=plot(x_fit,y_fit,'b','linewidth',1);
        xlim(XXLim);ylim(YYLim);
        title('All');
        text(0.95*XXLim(1)+0.05*XXLim(2),0.95*YYLim(2)+0.05*YYLim(1),{['P = ' num2str(stats.p(2),'%.2e')]});
        set(gca,'FontName','Arial','FontUnits','points','FontSize',6);
    end
end

fontsize(gcf, 7, 'points');
set(gcf,'FontName','Arial');
% exportgraphics(gcf,'figures/exFig10.eps','Colorspace','rgb','Resolution',600);
print(gcf, 'figures/exFig10.eps', '-depsc', '-vector');