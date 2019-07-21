<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            Taxonomy of sample with affy id ${affy_id}
        </h3>
        <table align="center" class="colored">
            %if len(taxonomies) == 0:
                <tr>
                    <td colspan="2">
                        This sample has no taxonomy
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Species Name</td>
                    <td>Genus Name</td>
                </tr>
                <% ctr = 0 %>
                %for taxonomy in taxonomies:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${taxonomy[0]}</td>
                        <td>${taxonomy[1]}</td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

