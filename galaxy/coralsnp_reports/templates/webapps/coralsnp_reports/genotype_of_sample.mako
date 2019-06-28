<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            Genotype of sample with affy id ${affy_id}
        </h3>
        <table align="center" class="colored">
            %if len(genotypes) == 0:
                <tr>
                    <td colspan="2">
                        This sample has no genotype
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Coral MLG Clonal ID</td>
                    <td>Coral MLG Rep Sample ID</td>
                    <td>Genetic Coral Species Call</td>
                </tr>
                <% ctr = 0 %>
                %for genotype in genotypes:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${genotype[0]}</td>
                        <td>${genotype[1]}</td>
                        <td>${genotype[2]}</td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

