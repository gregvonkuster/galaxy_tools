<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            Alleles of sample with affy id ${affy_id}
        </h3>
        <table align="center" class="colored">
            %if len(alleles) == 0:
                <tr>
                    <td colspan="2">
                        This sample has no alleles
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Allele</td>
                </tr>
                <% ctr = 0 %>
                %for allele in alleles:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${allele[1]}</td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

