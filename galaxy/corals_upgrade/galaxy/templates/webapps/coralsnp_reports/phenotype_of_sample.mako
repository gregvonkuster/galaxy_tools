<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            Phenotype of sample with affy id ${affy_id}
        </h3>
        <table align="center" class="colored">
            %if len(phenotypes) == 0:
                <tr>
                    <td colspan="2">
                        This sample has no phenotype
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Disease Resist</td>
                    <td>Bleach Resist</td>
                    <td>Mortality</td>
                    <td>TLE</td>
                    <td>Spawning</td>
                    <td>Sperm Motility</td>
                    <td>Healing Time</td>
                </tr>
                <% ctr = 0 %>
                %for phenotype in phenotypes:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${phenotype[0]}</td>
                        <td>${phenotype[1]}</td>
                        <td>${phenotype[2]}</td>
                        <td>${phenotype[3]}</td>
                        <td>${phenotype[4]}</td>
                        <td>${phenotype[5]}</td>
                        <td>${phenotype[6]}</td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

