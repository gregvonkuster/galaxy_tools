<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            Genotypes for sample with
            &nbsp;affy id&nbsp;${affy_id}&nbsp;
            &nbsp;sample id&nbsp;${sample_id}&nbsp;
            user specimen id&nbsp;${user_specimen_id}
        </h3>
        <table align="center" width="30%" class="colored">
            %if len(genotypes) == 0:
                <tr>
                    <td colspan="2">
                        There are no genotypes for this sample
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>coral mlg clonal id</td>
                    <td>coral mlg rep sample id</td>
                    <td>genetic coral species call</td>
                    <td>bcoral genet id</td>
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
                        <td>${genotype[3]}</td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>
